#!/usr/bin/env python3
"""Verify SHA-256 of exact original scientific/code files in the PACE repository."""
import csv,hashlib,sys
from pathlib import Path
root=Path(__file__).resolve().parent.parent
manifest=root/'MASTER_MANIFEST_SHA256.csv'
with manifest.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
errors=[]
for row in rows:
    p=root/row['relative_path']
    if not p.is_file(): errors.append('MISSING '+row['relative_path']);continue
    if p.stat().st_size != int(row['bytes']): errors.append('SIZE MISMATCH '+row['relative_path']);continue
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    if h.lower()!=row['sha256'].lower(): errors.append('HASH MISMATCH '+row['relative_path'])
print(f'Checked {len(rows)} scientific files, {len(errors)} errors')
for e in errors: print(e,file=sys.stderr)
if errors: raise SystemExit(1)
