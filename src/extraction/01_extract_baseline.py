"""Read-only baseline-only export; leave the frozen regression inputs intact."""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

import pandas as pd

ROOT = Path(os.environ['IR_REVISION_ROOT']).resolve()
PSQL = shutil.which('psql')
if not PSQL:
    raise RuntimeError('PostgreSQL psql client is required.')

for database, port, stay_key in [('mimiciii', 5433, 'icustay_id'), ('nwicu', 5432, 'stay_id')]:
    folder = ROOT / 'data/restricted'
    sql = ROOT / 'src/extraction' / database / '03_extract_baseline.sql'
    original = pd.read_csv(folder / f'{database}_patient_level.csv', usecols=['subject_id', 'hadm_id', stay_key])
    environment = os.environ.copy()
    environment['PGOPTIONS'] = '-c default_transaction_read_only=on -c work_mem=256MB'
    with tempfile.NamedTemporaryFile(dir=folder, suffix='.csv', delete=False) as out:
        staged = Path(out.name)
        os.chmod(staged, 0o600)
        subprocess.run([PSQL, '-X', '-w', '-v', 'ON_ERROR_STOP=1', '-h', 'localhost',
                        '-p', str(port), '-d', database, '--csv', '-f', str(sql)],
                       env=environment, stdout=out, check=True)
    data = pd.read_csv(staged)
    keys = ['subject_id', 'hadm_id', stay_key]
    assert not data[stay_key].duplicated().any()
    merged = original.merge(data, on=keys, validate='one_to_one', how='outer', indicator=True)
    assert merged['_merge'].eq('both').all(), 'Baseline export does not match the current eligible first-stay input.'
    staged.replace(folder / f'{database}_baseline_supplement.csv')
    print(f'{database}: {len(data):,} eligible stays; exact identifiers matched; baseline supplement saved.')
