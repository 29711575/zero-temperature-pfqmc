#!/usr/bin/env python3
import csv
import hashlib
from collections import defaultdict
from pathlib import Path

LV_ROOT = Path('/home/sunxr/new-pfqmc-zeroT-sign-20260827/reproduction/projector_kitaev/zeroT_sign_L6_L12_L18_20260827')
LT_ROOT = Path('/home/sunxr/new-pfqmc-zeroT-L-theta-20260827/reproduction/projector_kitaev/zeroT_sign_L_theta_20260827')
INPUTS = [('LV', LV_ROOT, LV_ROOT / 'per_seed.csv'),
          ('Ltheta', LT_ROOT, LT_ROOT / 'all_seeds_L_theta.csv')]
KEY_FIELDS = ('L', 'theta', 'V', 'seed', 'beta_trial', 'dt', 'delta', 'mu',
              'boundary', 'hs_scheme', 'burn', 'measurements',
              'provenance_executable_sha256')
OUT_FIELDS = (
    'task_index', 'source_dataset', 'archive_csv', 'archive_root',
    'archive_result_dir', 'archive_measurements', 'archive_measurements_sha256',
    'archive_source_commit', 'archive_executable_sha256', 'archive_job_id',
    'archive_status', 'archive_qc_status', 'L', 'theta', 'V', 'seed',
    'beta_trial', 'dt', 'delta', 'mu', 'boundary', 'hs_scheme', 'burn',
    'measurements', 'threads', 'diagnostic_stride', 'sign_stride',
    'legacy_average_sign', 'legacy_average_sign_err', 'legacy_runtime_seconds')


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    here = Path(__file__).resolve().parent
    rows = []
    for source, root, csv_path in INPUTS:
        with csv_path.open(newline='') as stream:
            for row in csv.DictReader(stream):
                row['_source'] = source
                row['_root'] = str(root)
                row['_csv'] = str(csv_path)
                rows.append(row)

    grouped = defaultdict(list)
    for row in rows:
        grouped[tuple(row.get(name, '') for name in KEY_FIELDS)].append(row)
    duplicates = sum(len(group) - 1 for group in grouped.values())
    if len(rows) != 186 or len(grouped) != 150 or duplicates != 36:
        raise RuntimeError(f'unexpected archive inventory: input={len(rows)} unique={len(grouped)} duplicates={duplicates}')

    output = []
    for group in grouped.values():
        # Prefer the original L/V row for the 36 theta=L references.
        row = sorted(group, key=lambda item: item['_source'] != 'LV')[0]
        root = Path(row['_root'])
        if row['_source'] == 'Ltheta' and row.get('source_dataset') == 'reused_theta_equals_L':
            root = LV_ROOT
        measurements = root / row['result_dir'] / 'measurements.csv'
        if not measurements.is_file():
            raise FileNotFoundError(measurements)
        if sum(1 for _ in measurements.open()) != int(row['measurements']) + 1:
            raise RuntimeError(f'wrong line count: {measurements}')
        output.append({
            'source_dataset': row['_source'],
            'archive_csv': row['_csv'],
            'archive_root': str(root),
            'archive_result_dir': row['result_dir'],
            'archive_measurements': str(measurements),
            'archive_measurements_sha256': sha256(measurements),
            'archive_source_commit': row['provenance_source_commit'],
            'archive_executable_sha256': row['provenance_executable_sha256'],
            'archive_job_id': row.get('provenance_job_id', ''),
            'archive_status': row.get('status', ''),
            'archive_qc_status': row.get('qc_status', ''),
            'L': row['L'], 'theta': row['theta'], 'V': row['V'],
            'seed': row['seed'], 'beta_trial': row['beta_trial'],
            'dt': row['dt'], 'delta': row['delta'], 'mu': row['mu'],
            'boundary': row['boundary'], 'hs_scheme': row['hs_scheme'],
            'burn': row['burn'], 'measurements': row['measurements'],
            'threads': row.get('threads', '1'),
            'diagnostic_stride': row.get('diagnostic_stride', '200'),
            'sign_stride': row.get('sign_recompute_stride', '20'),
            'legacy_average_sign': row['average_sign'],
            'legacy_average_sign_err': row['average_sign_err'],
            'legacy_runtime_seconds': row.get('runtime_seconds', ''),
        })
    output.sort(key=lambda r: (int(r['L']), float(r['theta']), int(r['boundary']),
                               int(r['hs_scheme']), float(r['V']), int(r['seed']),
                               r['archive_executable_sha256']))
    for index, row in enumerate(output):
        row['task_index'] = str(index)
    manifest = here / 'archive_manifest.csv'
    with manifest.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=OUT_FIELDS)
        writer.writeheader(); writer.writerows(output)
    counts = here / 'archive_inventory.txt'
    counts.write_text(
        f'input_rows={len(rows)}\nunique_seeds={len(output)}\n'
        f'duplicate_rows={duplicates}\npbc_hs0={sum(r["boundary"] == "0" and r["hs_scheme"] == "0" for r in output)}\n'
        f'obc_hs1={sum(r["boundary"] == "1" and r["hs_scheme"] == "1" for r in output)}\n')
    print(counts.read_text(), end='')


if __name__ == '__main__':
    main()
