#!/usr/bin/env python3
"""Verify published repository files (excluding this self-check and manifest)."""
import csv, hashlib
from pathlib import Path
root=Path(__file__).resolve().parent.parent
manifest=root/'PUBLIC_REPOSITORY_FILE_SHA256.csv'
errs=[];rows=[]
with manifest.open(encoding='utf-8',newline='') as fh:rows=list(csv.DictReader(fh))
for row in rows:
    p=root/row['relative_path']
    if not p.is_file():errs.append('MISSING '+row['relative_path']);continue
    actual=hashlib.sha256(p.read_bytes()).hexdigest()
    if actual!=row['sha256']:errs.append('MISMATCH '+row['relative_path'])
print(f'Checked {len(rows)} repository files; failures: {len(errs)}')
for e in errs:print(e)
if errs:raise SystemExit(1)
