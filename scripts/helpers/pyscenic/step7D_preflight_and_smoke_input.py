import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import loompy
from ctxcore.genesig import GeneSignature

expression_loom = Path(sys.argv[1])
regulon_dat = Path(sys.argv[2])
regulon_summary_file = Path(sys.argv[3])
smoke_loom = Path(sys.argv[4])
smoke_regulon_gmt = Path(sys.argv[5])
audit_file = Path(sys.argv[6])
expected_genes = int(sys.argv[7])
expected_cells = int(sys.argv[8])
expected_regulons = int(sys.argv[9])
smoke_cells_n = int(sys.argv[10])
smoke_regulons_n = int(sys.argv[11])
key_tfs = [x for x in sys.argv[12].split(',') if x]

with regulon_dat.open('rb') as handle:
    regulons = pickle.load(handle)

if len(regulons) != expected_regulons:
    raise RuntimeError(f'Expected {expected_regulons} regulons, observed {len(regulons)}')

regulon_names = [str(x.name) for x in regulons]
if len(set(regulon_names)) != len(regulon_names):
    raise RuntimeError('Duplicated regulon names detected')

regulon_tfs = [
    str(getattr(x, 'transcription_factor', str(x.name).split('(')[0]))
    for x in regulons
]

regulon_sizes = [len(x.genes) for x in regulons]
regulon_summary = pd.read_csv(regulon_summary_file, sep='\t')

with loompy.connect(str(expression_loom), mode='r') as ds:
    loom_shape = tuple(ds.shape)
    genes = [str(x) for x in ds.ra.Gene]
    cells = [str(x) for x in ds.ca.CellID]
    branch_states = [str(x) for x in ds.ca.BranchState]

    if loom_shape != (expected_genes, expected_cells):
        raise RuntimeError(f'Unexpected expression loom shape: {loom_shape}')

    if len(set(genes)) != len(genes):
        raise RuntimeError('Duplicated genes in expression loom')

    if len(set(cells)) != len(cells):
        raise RuntimeError('Duplicated cells in expression loom')

    expected_states = {'Shared_prebranch': 116, 'Arm_1': 2289, 'Arm_2': 1798}
    observed_states = {state: branch_states.count(state) for state in expected_states}

    if observed_states != expected_states:
        raise RuntimeError(f'Frozen state counts differ: {observed_states}')

    loom_gene_set = set(genes)
    missing_targets = sorted(
        {gene for regulon in regulons for gene in regulon.genes} - loom_gene_set
    )

    if missing_targets:
        raise RuntimeError(f'{len(missing_targets)} regulon targets missing from loom')

    state_arrays = {state: np.flatnonzero(np.asarray(branch_states) == state) for state in expected_states}
    selected_cells = []
    per_state = max(1, smoke_cells_n // 3)

    for state in ['Shared_prebranch', 'Arm_1', 'Arm_2']:
        indices = state_arrays[state]
        selected_cells.extend(indices[:min(per_state, len(indices))].tolist())

    selected_cells = sorted(set(selected_cells))
    if len(selected_cells) < smoke_cells_n:
        remaining = [i for i in range(len(cells)) if i not in set(selected_cells)]
        selected_cells.extend(remaining[:smoke_cells_n - len(selected_cells)])

    selected_cells = np.asarray(sorted(selected_cells[:smoke_cells_n]), dtype=int)
    smoke_matrix = ds[:, selected_cells]
    smoke_cell_ids = np.asarray(ds.ca.CellID[selected_cells], dtype=object)

if smoke_loom.exists():
    smoke_loom.unlink()

loompy.create(
    str(smoke_loom),
    smoke_matrix,
    {'Gene': np.asarray(genes, dtype=object)},
    {'CellID': smoke_cell_ids}
)

preferred = []
for tf in key_tfs:
    for index, regulon_tf in enumerate(regulon_tfs):
        if regulon_tf == tf and index not in preferred:
            preferred.append(index)

size_order = list(np.argsort(np.asarray(regulon_sizes))[::-1])
selected_regulon_indices = []

for index in preferred + size_order:
    if index not in selected_regulon_indices:
        selected_regulon_indices.append(index)
    if len(selected_regulon_indices) >= smoke_regulons_n:
        break

smoke_regulons = [regulons[i] for i in selected_regulon_indices]

GeneSignature.to_gmt(str(smoke_regulon_gmt), smoke_regulons)

primary_regulons = sum(size >= 10 for size in regulon_sizes)
small_regulons = sum(size < 10 for size in regulon_sizes)
key_tfs_present = len(set(regulon_tfs).intersection(key_tfs))

metrics = [
    ('loom_genes', len(genes)),
    ('loom_cells', len(cells)),
    ('state_Shared_prebranch', observed_states['Shared_prebranch']),
    ('state_Arm_1', observed_states['Arm_1']),
    ('state_Arm_2', observed_states['Arm_2']),
    ('regulon_count', len(regulons)),
    ('unique_regulon_TFs', len(set(regulon_tfs))),
    ('primary_regulons_n_targets_ge_10', primary_regulons),
    ('small_regulons_n_targets_lt_10', small_regulons),
    ('minimum_regulon_targets', min(regulon_sizes)),
    ('median_regulon_targets', float(np.median(regulon_sizes))),
    ('maximum_regulon_targets', max(regulon_sizes)),
    ('regulon_targets_missing_from_loom', len(missing_targets)),
    ('key_TFs_with_regulons', key_tfs_present),
    ('smoke_cells', len(selected_cells)),
    ('smoke_regulons', len(smoke_regulons))
]

with audit_file.open('w', encoding='utf-8', newline='\n') as handle:
    handle.write('metric\tvalue\n')
    for metric, value in metrics:
        handle.write(f'{metric}\t{value}\n')

print('STEP7D_PREFLIGHT_PASS')
for metric, value in metrics:
    print(f'{metric}: {value}')
