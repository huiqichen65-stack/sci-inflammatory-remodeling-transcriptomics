import gzip
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import loompy

input_loom = Path(sys.argv[1])
aucell_loom = Path(sys.argv[2])
regulon_dat = Path(sys.argv[3])
auc_tsv = Path(sys.argv[4])
matrix_audit_file = Path(sys.argv[5])
distribution_file = Path(sys.argv[6])
state_summary_file = Path(sys.argv[7])
key_tf_file = Path(sys.argv[8])
timepoint_summary_file = Path(sys.argv[9])
expected_cells = int(sys.argv[10])
expected_regulons = int(sys.argv[11])
primary_min_targets = int(sys.argv[12])
min_unique_values = int(sys.argv[13])
low_sd = float(sys.argv[14])
low_range = float(sys.argv[15])
key_tfs = [x for x in sys.argv[16].split(',') if x]

with regulon_dat.open('rb') as handle:
    regulons = pickle.load(handle)

regulon_names_expected = [str(x.name) for x in regulons]
regulon_tfs = [
    str(getattr(x, 'transcription_factor', str(x.name).split('(')[0]))
    for x in regulons
]
regulon_sizes = [len(x.genes) for x in regulons]

with loompy.connect(str(input_loom), mode='r') as source:
    source_cells = [str(x) for x in source.ca.CellID]
    metadata = pd.DataFrame({
        'CellID': source_cells,
        'BranchState': [str(x) for x in source.ca.BranchState],
        'SampleID': [str(x) for x in source.ca.SampleID],
        'Timepoint': [str(x) for x in source.ca.Timepoint],
        'Pseudotime': np.asarray(source.ca.Pseudotime, dtype=float),
        'Y29SignedProgress': np.asarray(source.ca.Y29SignedProgress, dtype=float)
    })

with loompy.connect(str(aucell_loom), mode='r', validate=False) as ds:
    output_cells = [str(x) for x in ds.ca.CellID]

    auc_attribute = None
    for candidate in ['RegulonsAUC', 'RegulonAUC']:
        if candidate in ds.ca:
            auc_attribute = candidate
            break

    if auc_attribute is None:
        raise RuntimeError(
            'No RegulonsAUC/RegulonAUC column attribute found. Available: ' +
            ', '.join(ds.ca.keys())
        )

    raw_auc = np.asarray(ds.ca[auc_attribute])

if output_cells != source_cells:
    raise RuntimeError('AUCell loom cell order differs from frozen input loom')

if raw_auc.dtype.names is None:
    raise RuntimeError(
        f'AUC attribute is not a named structured array: shape={raw_auc.shape}, dtype={raw_auc.dtype}'
    )

auc_dict = {}
raw_to_normalized_names = {}

for raw_name in raw_auc.dtype.names:
    raw_name_string = str(raw_name)
    normalized_name = raw_name_string[:-1] if raw_name_string.endswith(',') else raw_name_string

    if normalized_name in auc_dict:
        raise RuntimeError(
            f'Regulon-name normalization collision: {raw_name_string} -> {normalized_name}'
        )

    values = np.asarray(raw_auc[raw_name]).reshape(-1).astype(float)
    auc_dict[normalized_name] = values
    raw_to_normalized_names[raw_name_string] = normalized_name

auc = pd.DataFrame(auc_dict, index=output_cells)
auc.index.name = 'CellID'

if auc.shape != (expected_cells, expected_regulons):
    raise RuntimeError(f'Unexpected AUC matrix shape: {auc.shape}')

if set(auc.columns) != set(regulon_names_expected):
    missing = sorted(set(regulon_names_expected) - set(auc.columns))
    extra = sorted(set(auc.columns) - set(regulon_names_expected))
    raise RuntimeError(f'Regulon-name mismatch. Missing={missing[:10]}, extra={extra[:10]}')

auc = auc.loc[:, regulon_names_expected]

values = auc.to_numpy(dtype=float)
finite_mask = np.isfinite(values)
bad_values = int((~finite_mask).sum())
negative_values = int((values < -1e-12).sum())
above_one_values = int((values > 1.000001).sum())

if bad_values > 0:
    raise RuntimeError(f'AUC matrix contains {bad_values} non-finite values')

auc.reset_index().to_csv(
    auc_tsv,
    sep='\t',
    index=False,
    compression='gzip'
)

distribution_rows = []

for index, regulon_name in enumerate(regulon_names_expected):
    vector = auc[regulon_name].to_numpy(dtype=float)
    quantiles = np.quantile(vector, [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99])
    unique_values = int(np.unique(vector).size)
    standard_deviation = float(np.std(vector, ddof=1))
    observed_range = float(np.max(vector) - np.min(vector))
    low_information = (
        standard_deviation <= low_sd
        or observed_range <= low_range
        or unique_values < min_unique_values
    )

    distribution_rows.append({
        'regulon': regulon_name,
        'TF': regulon_tfs[index],
        'n_targets': regulon_sizes[index],
        'primary_n_targets_ge_10': regulon_sizes[index] >= primary_min_targets,
        'mean_AUC': float(np.mean(vector)),
        'sd_AUC': standard_deviation,
        'minimum_AUC': float(np.min(vector)),
        'q01_AUC': float(quantiles[0]),
        'q05_AUC': float(quantiles[1]),
        'q25_AUC': float(quantiles[2]),
        'median_AUC': float(quantiles[3]),
        'q75_AUC': float(quantiles[4]),
        'q95_AUC': float(quantiles[5]),
        'q99_AUC': float(quantiles[6]),
        'maximum_AUC': float(np.max(vector)),
        'IQR_AUC': float(quantiles[4] - quantiles[2]),
        'range_AUC': observed_range,
        'zero_fraction': float(np.mean(vector == 0)),
        'unique_AUC_values': unique_values,
        'low_information_flag': bool(low_information)
    })

distribution = pd.DataFrame(distribution_rows)
distribution = distribution.sort_values(
    ['low_information_flag', 'sd_AUC', 'n_targets'],
    ascending=[True, False, False]
)
distribution.to_csv(distribution_file, sep='\t', index=False)

auc_with_metadata = metadata.set_index('CellID').join(auc, how='left')

state_rows = []
for state, frame in auc_with_metadata.groupby('BranchState', sort=False):
    for regulon_name in regulon_names_expected:
        vector = frame[regulon_name].to_numpy(dtype=float)
        state_rows.append({
            'BranchState': state,
            'regulon': regulon_name,
            'TF': regulon_tfs[regulon_names_expected.index(regulon_name)],
            'n_cells': int(frame.shape[0]),
            'mean_AUC': float(np.mean(vector)),
            'median_AUC': float(np.median(vector)),
            'sd_AUC': float(np.std(vector, ddof=1)) if len(vector) > 1 else np.nan
        })

pd.DataFrame(state_rows).to_csv(state_summary_file, sep='\t', index=False)

time_rows = []
for timepoint, frame in auc_with_metadata.groupby('Timepoint', sort=False):
    for regulon_name in regulon_names_expected:
        vector = frame[regulon_name].to_numpy(dtype=float)
        time_rows.append({
            'Timepoint': timepoint,
            'regulon': regulon_name,
            'TF': regulon_tfs[regulon_names_expected.index(regulon_name)],
            'n_cells': int(frame.shape[0]),
            'mean_AUC': float(np.mean(vector)),
            'median_AUC': float(np.median(vector)),
            'sd_AUC': float(np.std(vector, ddof=1)) if len(vector) > 1 else np.nan
        })

pd.DataFrame(time_rows).to_csv(timepoint_summary_file, sep='\t', index=False)

distribution_indexed = distribution.set_index('TF', drop=False)
key_rows = []

for tf in key_tfs:
    matching = [i for i, value in enumerate(regulon_tfs) if value == tf]

    if not matching:
        key_rows.append({
            'TF': tf,
            'regulon': '',
            'present_in_regulons': False,
            'n_targets': 0,
            'primary_n_targets_ge_10': False,
            'mean_AUC': np.nan,
            'sd_AUC': np.nan,
            'median_AUC': np.nan,
            'q95_AUC': np.nan,
            'maximum_AUC': np.nan,
            'unique_AUC_values': 0,
            'low_information_flag': np.nan
        })
        continue

    for index in matching:
        regulon_name = regulon_names_expected[index]
        row = distribution.loc[distribution['regulon'] == regulon_name].iloc[0]
        key_rows.append({
            'TF': tf,
            'regulon': regulon_name,
            'present_in_regulons': True,
            'n_targets': int(row['n_targets']),
            'primary_n_targets_ge_10': bool(row['primary_n_targets_ge_10']),
            'mean_AUC': float(row['mean_AUC']),
            'sd_AUC': float(row['sd_AUC']),
            'median_AUC': float(row['median_AUC']),
            'q95_AUC': float(row['q95_AUC']),
            'maximum_AUC': float(row['maximum_AUC']),
            'unique_AUC_values': int(row['unique_AUC_values']),
            'low_information_flag': bool(row['low_information_flag'])
        })

pd.DataFrame(key_rows).to_csv(key_tf_file, sep='\t', index=False)

low_information_count = int(distribution['low_information_flag'].sum())
primary_count = int(distribution['primary_n_targets_ge_10'].sum())
nonzero_regulons = int((distribution['maximum_AUC'] > 0).sum())

metrics = [
    ('AUC_cells', int(auc.shape[0])),
    ('AUC_regulons', int(auc.shape[1])),
    ('AUC_values', int(values.size)),
    ('non_finite_AUC_values', bad_values),
    ('negative_AUC_values', negative_values),
    ('AUC_values_above_1', above_one_values),
    ('minimum_AUC', float(values.min())),
    ('median_AUC_all_values', float(np.median(values))),
    ('mean_AUC_all_values', float(values.mean())),
    ('maximum_AUC', float(values.max())),
    ('primary_regulons_n_targets_ge_10', primary_count),
    ('small_regulons_n_targets_lt_10', int(expected_regulons - primary_count)),
    ('nonzero_regulons', nonzero_regulons),
    ('low_information_regulons', low_information_count),
    ('informative_regulons', int(expected_regulons - low_information_count)),
    ('low_information_fraction', low_information_count / expected_regulons),
    ('cell_order_exact_match', True),
    ('regulon_name_set_exact_match', True),
    ('raw_regulon_names_with_terminal_comma', int(sum(name.endswith(',') for name in raw_to_normalized_names))),
    ('regulon_name_normalization_collisions', 0),
    ('AUC_tsv_size_bytes', auc_tsv.stat().st_size)
]

pd.DataFrame(metrics, columns=['metric', 'value']).to_csv(
    matrix_audit_file,
    index=False
)

print('STEP7D_AUC_EXTRACTION_PASS')
for metric, value in metrics:
    print(f'{metric}: {value}')
