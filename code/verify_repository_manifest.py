#!/usr/bin/env python3
"""Verify public source and documentation files other than the manifests themselves."""
import csv,hashlib,sys
from pathlib import Path
root=Path(__file__).resolve().parent.parent
manifest=root/'PUBLIC_REPOSITORY_FILE_SHA256.csv'
with manifest.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
errors=[]
for row in rows:
    p=root/row['relative_path']
    if not p.is_file(): errors.append('MISSING '+row['relative_path']);continue
    if p.stat().st_size!=int(row['bytes']): errors.append('SIZE MISMATCH '+row['relative_path']);continue
    if hashlib.sha256(p.read_bytes()).hexdigest().lower()!=row['sha256'].lower(): errors.append('HASH MISMATCH '+row['relative_path'])
print(f'Checked {len(rows)} public repository files, {len(errors)} errors')
for e in errors: print(e,file=sys.stderr)
if errors: raise SystemExit(1)
