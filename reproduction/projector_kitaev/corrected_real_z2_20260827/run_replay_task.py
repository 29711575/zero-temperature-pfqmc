#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('manifest')
    parser.add_argument('task_index', type=int)
    parser.add_argument('output_root')
    parser.add_argument('--mp-spot-stride', type=int, default=0)
    args = parser.parse_args()
    base = Path(__file__).resolve().parent
    with open(args.manifest, newline='') as stream:
        rows = list(csv.DictReader(stream))
    row = rows[args.task_index]
    if int(row['task_index']) != args.task_index:
        raise RuntimeError('manifest task index mismatch')
    out = Path(args.output_root) / f'task_{args.task_index:03d}'
    out.mkdir(parents=True, exist_ok=False)
    executable = base.parent / 'real_z2_raw_checker_fix_20260827/bin/projector_real_z2_driver'
    command = [str(executable), row['L'], row['theta'], row['beta_trial'], row['dt'],
               row['V'], row['delta'], row['mu'], row['boundary'], row['hs_scheme'],
               row['seed'], row['burn'], row['measurements'], '1',
               str(out / 'comparison.csv'), row['diagnostic_stride'], row['sign_stride'],
               str(out / 'legacy_replay.csv'), str(args.mp_spot_stride)]
    environment = dict(os.environ)
    environment.update(OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       OPENBLAS_NUM_THREADS='1')
    started = time.monotonic()
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=environment)
    elapsed = time.monotonic() - started
    (out / 'stdout.log').write_text(completed.stdout)
    (out / 'stderr.log').write_text(completed.stderr)
    (out / 'command.txt').write_text(' '.join(command) + '\n')
    result = None
    for line in reversed(completed.stdout.splitlines()):
        try:
            result = json.loads(line)
            break
        except json.JSONDecodeError:
            pass
    legacy = out / 'legacy_replay.csv'
    legacy_hash = sha256(legacy) if legacy.is_file() else None
    expected_hash = row['archive_measurements_sha256']
    wrapper = {
        'task_index': args.task_index,
        'exit_code': completed.returncode,
        'wrapper_elapsed_seconds': elapsed,
        'archive_measurements_sha256': expected_hash,
        'legacy_replay_sha256': legacy_hash,
        'archive_trajectory_exact': legacy_hash == expected_hash,
        'new_executable_sha256': sha256(executable),
        'mp_spot_stride': args.mp_spot_stride,
        'result_json_present': result is not None,
        'result': result,
    }
    (out / 'wrapper.json').write_text(json.dumps(wrapper, indent=2, sort_keys=True) + '\n')
    # Keep the PBS array healthy so every archived failure is available to classification.
    print(json.dumps({key: wrapper[key] for key in ('task_index', 'exit_code',
          'archive_trajectory_exact', 'result_json_present')}))


if __name__ == '__main__':
    main()
