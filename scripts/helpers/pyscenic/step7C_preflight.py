import sys
from pathlib import Path

import numpy as np
import pandas as pd
import loompy
from ctxcore.rnkdb import FeatherRankingDatabase

adjacency_file = Path(sys.argv[1])
loom_file = Path(sys.argv[2])
tf_file = Path(sys.argv[3])
annotation_file = Path(sys.argv[4])
feather_file = Path(sys.argv[5])
audit_file = Path(sys.argv[6])
expected_genes = int(sys.argv[7])
expected_cells = int(sys.argv[8])

adj = pd.read_csv(adjacency_file, sep='\t', dtype={'TF': str, 'target': str})
required = {'TF', 'target', 'importance'}
if not required.issubset(adj.columns):
    raise RuntimeError(f'Adjacency columns missing: {required - set(adj.columns)}')

with loompy.connect(str(loom_file), mode='r') as ds:
    loom_shape = tuple(ds.shape)
    genes = [str(x) for x in ds.ra.Gene]
    cells = [str(x) for x in ds.ca.CellID]
    branch_states = [str(x) for x in ds.ca.BranchState]

if loom_shape != (expected_genes, expected_cells):
    raise RuntimeError(f'Unexpected loom shape: {loom_shape}')
if len(set(genes)) != len(genes):
    raise RuntimeError('Duplicated loom genes detected')
if len(set(cells)) != len(cells):
    raise RuntimeError('Duplicated loom cells detected')

expected_states = {'Shared_prebranch': 116, 'Arm_1': 2289, 'Arm_2': 1798}
observed_states = {state: branch_states.count(state) for state in expected_states}
if observed_states != expected_states:
    raise RuntimeError(f'Frozen branch counts differ: {observed_states}')

network_genes = set(adj['TF']).union(set(adj['target']))
missing_network_genes = sorted(network_genes - set(genes))
if missing_network_genes:
    raise RuntimeError(f'{len(missing_network_genes)} network genes missing from loom')

with tf_file.open('r', encoding='utf-8-sig') as handle:
    tf_list = [line.strip() for line in handle if line.strip()]
tf_overlap = sorted(set(adj['TF']).intersection(tf_list))

database = FeatherRankingDatabase(fname=str(feather_file), name='mm10_10kb')
db_genes = set(database.genes)
db_overlap = set(genes).intersection(db_genes)

annotation = pd.read_csv(annotation_file, sep='\t', low_memory=False)
if annotation.shape[0] < 1000:
    raise RuntimeError(f'Motif annotation unexpectedly small: {annotation.shape}')

metrics = [
    ('loom_genes', len(genes)),
    ('loom_cells', len(cells)),
    ('state_Shared_prebranch', observed_states['Shared_prebranch']),
    ('state_Arm_1', observed_states['Arm_1']),
    ('state_Arm_2', observed_states['Arm_2']),
    ('adjacency_edges', int(adj.shape[0])),
    ('adjacency_TFs', int(adj['TF'].nunique())),
    ('adjacency_targets', int(adj['target'].nunique())),
    ('network_genes_missing_from_loom', len(missing_network_genes)),
    ('GRN_TFs_in_mouse_TF_list', len(tf_overlap)),
    ('feather_database_genes', len(db_genes)),
    ('loom_feather_overlap_genes', len(db_overlap)),
    ('loom_feather_overlap_fraction', len(db_overlap) / len(genes)),
    ('motif_annotation_rows', int(annotation.shape[0])),
    ('motif_annotation_columns', int(annotation.shape[1]))
]

with audit_file.open('w', encoding='utf-8', newline='\n') as handle:
    handle.write('metric\tvalue\n')
    for metric, value in metrics:
        handle.write(f'{metric}\t{value}\n')

print('STEP7C_PREFLIGHT_PASS')
for metric, value in metrics:
    print(f'{metric}: {value}')
