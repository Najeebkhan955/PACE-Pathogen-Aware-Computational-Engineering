#!/usr/bin/env python3
"""Check unchanged scientific files listed in MASTER_MANIFEST_SHA256.csv."""
import csv, hashlib
from pathlib import Path
root=Path(__file__).resolve().parent.parent
manifest=root/'MASTER_MANIFEST_SHA256.csv'
errors=[]
with manifest.open(encoding='utf-8',newline='') as f:
    rows=list(csv.DictReader(f))
for row in rows:
    path=root/row['relative_path']
    if not path.is_file(): errors.append('MISSING '+row['relative_path']); continue
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    if actual!=row['SHA256']: errors.append('SHA256 MISMATCH '+row['relative_path'])
print(f'Checked {len(rows)} scientific files; failures: {len(errors)}')
for e in errors: print(e)
if errors:raise SystemExit(1)
