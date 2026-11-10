import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import loompy

from ctxcore.genesig import GeneSignature
from pyscenic.transform import df2regulons
from pyscenic.utils import load_motifs

motif_file = Path(sys.argv[1])
adjacency_file = Path(sys.argv[2])
loom_file = Path(sys.argv[3])
dat_file = Path(sys.argv[4])
gmt_file = Path(sys.argv[5])
json_file = Path(sys.argv[6])
metrics_file = Path(sys.argv[7])
regulon_summary_file = Path(sys.argv[8])
targets_file = Path(sys.argv[9])
key_tf_file = Path(sys.argv[10])
motif_tf_file = Path(sys.argv[11])
key_tfs = [x for x in sys.argv[12].split(',') if x]

motifs = load_motifs(str(motif_file), sep=',')
if motifs.empty:
    raise RuntimeError('Motif enrichment table is empty')

regulons = list(df2regulons(motifs))
if len(regulons) == 0:
    raise RuntimeError('No regulons were derived from the motif table')

with loompy.connect(str(loom_file), mode='r') as ds:
    loom_genes = set(str(x) for x in ds.ra.Gene)

adj = pd.read_csv(
    adjacency_file,
    sep='\t',
    dtype={'TF': str, 'target': str}
)

grn_summary = (
    adj.groupby('TF', sort=False)
       .agg(
           GRN_n_targets=('target', 'nunique'),
           GRN_total_importance=('importance', 'sum')
       )
       .reset_index()
)

summary_rows = []
target_rows = []
name_to_targets = {}

for regulon in regulons:
    name = str(regulon.name)
    tf = str(getattr(regulon, 'transcription_factor', name.split('(')[0]))
    gene2weight = dict(regulon.gene2weight)
    targets = list(gene2weight.keys())
    weights = [float(gene2weight[g]) for g in targets]
    context = sorted(str(x) for x in getattr(regulon, 'context', []))
    sign = 'negative' if name.endswith('(-)') else 'positive'

    name_to_targets[name] = targets

    summary_rows.append({
        'regulon': name,
        'TF': tf,
        'sign': sign,
        'n_targets': len(targets),
        'total_weight': float(np.sum(weights)) if weights else 0.0,
        'mean_weight': float(np.mean(weights)) if weights else np.nan,
        'median_weight': float(np.median(weights)) if weights else np.nan,
        'maximum_weight': float(np.max(weights)) if weights else np.nan,
        'context': ' | '.join(context)
    })

    for target in targets:
        target_rows.append({
            'regulon': name,
            'TF': tf,
            'sign': sign,
            'target': str(target),
            'weight': float(gene2weight[target])
        })

regulon_summary = pd.DataFrame(summary_rows)
targets_long = pd.DataFrame(target_rows)

if regulon_summary['regulon'].duplicated().any():
    raise RuntimeError('Duplicated regulon names detected')

targets_outside_loom = sorted(set(targets_long['target']) - loom_genes)
duplicate_regulon_targets = int(
    targets_long.duplicated(['regulon', 'target']).sum()
)

with dat_file.open('wb') as handle:
    pickle.dump(regulons, handle)

GeneSignature.to_gmt(str(gmt_file), regulons)

with json_file.open('w', encoding='utf-8') as handle:
    json.dump(name_to_targets, handle, ensure_ascii=False, indent=2)

regulon_summary = regulon_summary.sort_values(
    ['n_targets', 'total_weight'],
    ascending=[False, False]
)
regulon_summary.to_csv(regulon_summary_file, sep='\t', index=False)
targets_long.to_csv(targets_file, sep='\t', index=False)

regulon_tf_summary = (
    regulon_summary.groupby('TF', sort=False)
    .agg(
        regulon_count=('regulon', 'nunique'),
        pruned_target_count=('n_targets', 'sum'),
        maximum_single_regulon_targets=('n_targets', 'max'),
        regulon_names=('regulon', lambda x: ' | '.join(sorted(set(x))))
    )
    .reset_index()
)

key_table = pd.DataFrame({'TF': key_tfs})
key_table = key_table.merge(grn_summary, on='TF', how='left')
key_table = key_table.merge(regulon_tf_summary, on='TF', how='left')
key_table['present_in_GRN'] = key_table['GRN_n_targets'].notna()
key_table['present_in_pruned_regulon'] = key_table['regulon_count'].notna()
key_table['pruned_target_count'] = key_table['pruned_target_count'].fillna(0).astype(int)
key_table['regulon_count'] = key_table['regulon_count'].fillna(0).astype(int)
key_table['retention_fraction_vs_GRN_targets'] = np.where(
    key_table['GRN_n_targets'].fillna(0) > 0,
    key_table['pruned_target_count'] / key_table['GRN_n_targets'],
    np.nan
)
key_table.to_csv(key_tf_file, sep='\t', index=False)

motif_reset = motifs.reset_index()
flattened_columns = []
for col in motif_reset.columns:
    if isinstance(col, tuple):
        parts = [str(x) for x in col if str(x) and not str(x).startswith('Unnamed:')]
        flattened_columns.append('__'.join(parts))
    else:
        flattened_columns.append(str(col))
motif_reset.columns = flattened_columns

tf_column = next(
    (c for c in motif_reset.columns if c == 'TF' or c.endswith('__TF')),
    motif_reset.columns[0]
)
nes_column = next(
    (c for c in motif_reset.columns if c == 'NES' or c.endswith('__NES')),
    None
)

if nes_column is not None:
    motif_tf_summary = (
        motif_reset.groupby(tf_column)
        .agg(
            enriched_motif_rows=(tf_column, 'size'),
            maximum_NES=(nes_column, 'max'),
            median_NES=(nes_column, 'median')
        )
        .reset_index()
        .rename(columns={tf_column: 'TF'})
    )
else:
    motif_tf_summary = (
        motif_reset.groupby(tf_column)
        .size()
        .reset_index(name='enriched_motif_rows')
        .rename(columns={tf_column: 'TF'})
    )

motif_tf_summary = motif_tf_summary.sort_values(
    'enriched_motif_rows',
    ascending=False
)
motif_tf_summary.to_csv(motif_tf_file, sep='\t', index=False)

unique_regulon_tfs = int(regulon_summary['TF'].nunique())
unique_pruned_targets = int(targets_long['target'].nunique())
positive_regulons = int((regulon_summary['sign'] == 'positive').sum())
negative_regulons = int((regulon_summary['sign'] == 'negative').sum())
key_tfs_retained = int(key_table['present_in_pruned_regulon'].sum())

metrics = [
    ('motif_rows', int(motifs.shape[0])),
    ('motif_table_unique_TFs', int(motifs.index.get_level_values(0).nunique())),
    ('regulon_count', int(regulon_summary.shape[0])),
    ('unique_regulon_TFs', unique_regulon_tfs),
    ('positive_regulons', positive_regulons),
    ('negative_regulons', negative_regulons),
    ('pruned_regulon_target_edges', int(targets_long.shape[0])),
    ('unique_pruned_targets', unique_pruned_targets),
    ('minimum_regulon_targets', int(regulon_summary['n_targets'].min())),
    ('median_regulon_targets', float(regulon_summary['n_targets'].median())),
    ('mean_regulon_targets', float(regulon_summary['n_targets'].mean())),
    ('maximum_regulon_targets', int(regulon_summary['n_targets'].max())),
    ('original_GRN_edges', int(adj.shape[0])),
    ('original_GRN_TFs', int(adj['TF'].nunique())),
    ('TF_retention_fraction', unique_regulon_tfs / adj['TF'].nunique()),
    ('edge_retention_fraction', targets_long.shape[0] / adj.shape[0]),
    ('targets_outside_loom', len(targets_outside_loom)),
    ('duplicate_regulon_target_pairs', duplicate_regulon_targets),
    ('key_TFs_retained', key_tfs_retained),
    ('dat_file_size_bytes', dat_file.stat().st_size),
    ('gmt_file_size_bytes', gmt_file.stat().st_size),
    ('json_file_size_bytes', json_file.stat().st_size)
]

with metrics_file.open('w', encoding='utf-8', newline='\n') as handle:
    handle.write('metric\tvalue\n')
    for metric, value in metrics:
        handle.write(f'{metric}\t{value}\n')

print('STEP7C_REGULON_AUDIT_PASS')
for metric, value in metrics:
    print(f'{metric}: {value}')
