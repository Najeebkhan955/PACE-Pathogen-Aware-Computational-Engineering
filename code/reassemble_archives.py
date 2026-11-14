#!/usr/bin/env python3
"""Rebuild and verify losslessly split S06/S10 research archives.

Usage: python code/reassemble_archives.py
Optional: python code/reassemble_archives.py --verify-only

Requires only Python 3 standard library. Never modifies the Git-tracked part files.
"""
import argparse
import csv
import hashlib
from pathlib import Path

CHUNK_SIZE = 1024 * 1024


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(CHUNK_SIZE), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-only', action='store_true', help='Verify uploaded parts without writing reconstructed archives')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    manifest = root / 'LARGE_ARCHIVE_PARTS_SHA256.csv'
    with manifest.open(newline='', encoding='utf-8-sig') as stream:
        records = list(csv.DictReader(stream))
    archives = {}
    for row in records:
        archives.setdefault((row['group'], row['archive']), []).append(row)
    if not archives:
        raise RuntimeError('No archive parts listed in manifest')
    for (group, name), rows in sorted(archives.items()):
        rows.sort(key=lambda r: int(r['part_order']))
        if len(rows) != int(rows[0]['total_parts']):
            raise RuntimeError(f'{name}: not all expected parts are listed')
        if [int(r['part_order']) for r in rows] != list(range(1, len(rows) + 1)):
            raise RuntimeError(f'{name}: nonconsecutive parts')
        expected = rows[0]['reconstructed_archive_sha256']
        target = root / group / 'DATA' / name
        if target.exists():
            observed = digest(target)
            if observed != expected:
                raise RuntimeError(f'{target}: existing archive hash mismatch, refusing to overwrite')
            print(f'OK existing archive: {target.relative_to(root)}; SHA256 {observed}')
            continue
        whole = hashlib.sha256()
        total = 0
        if not args.verify_only:
            target.parent.mkdir(parents=True, exist_ok=True)
            temp = target.with_name(target.name + '.building')
            output = temp.open('wb')
        else:
            output = None
        try:
            for row in rows:
                path = root / row['part_path']
                if not path.is_file():
                    raise FileNotFoundError(f'Missing part: {row["part_path"]}')
                if path.stat().st_size != int(row['part_bytes']) or digest(path) != row['part_sha256']:
                    raise RuntimeError(f'Part failed integrity check: {row["part_path"]}')
                with path.open('rb') as stream:
                    for block in iter(lambda: stream.read(CHUNK_SIZE), b''):
                        whole.update(block)
                        total += len(block)
                        if output:
                            output.write(block)
        except Exception:
            if output:
                output.close()
                temp.unlink(missing_ok=True)
            raise
        if output:
            output.close()
        observed = whole.hexdigest()
        if total != int(rows[0]['archive_bytes']) or observed != expected:
            if output:
                temp.unlink(missing_ok=True)
            raise RuntimeError(f'{name}: reconstructed archive hash/size mismatch')
        if output:
            temp.replace(target)
        print(f'OK {"verified parts" if args.verify_only else "reconstructed"}: {group}/DATA/{name}; {total:,} bytes; SHA256 {observed}')
    print(f'SUCCESS: {len(archives)} research archives validated.')


if __name__ == '__main__':
    main()
