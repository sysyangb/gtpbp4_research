#!/usr/bin/env python3
"""Restore retained files to their original relative layout in a new directory.

Copies (never links) files so running old analysis scripts cannot alter this archive.
An optional original-root prefix can be replaced in the copied code. Files in
this repository remain unchanged.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path,
                        help='New, empty directory for the restored analysis project')
    parser.add_argument('--original-root-prefix',
                        help='Optional absolute path prefix to replace in copied code')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    dest = args.output.expanduser().resolve()
    if dest == root or root in dest.parents:
        parser.error('Use a directory outside the archive to protect saved results.')
    if dest.exists() and (not dest.is_dir() or any(dest.iterdir())):
        parser.error('Output must be absent or empty; existing files will not be overwritten.')
    metadata = json.loads((root / 'provenance/archive_metadata.json').read_text())
    original_prefix = args.original_root_prefix or metadata.get('original_project_root')
    if original_prefix is not None:
        if not isinstance(original_prefix, str) or not original_prefix.strip():
            parser.error('Original-root prefix must be a nonempty path when supplied.')
        if not Path(original_prefix).is_absolute():
            parser.error('Original-root prefix must be an absolute path.')
    with (root / 'provenance/file_manifest.tsv').open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    verified = set()
    for row in rows:
        for field in ('source_path', 'archive_path'):
            rel = Path(row[field])
            if not row[field] or rel.is_absolute() or '..' in rel.parts:
                raise SystemExit('Unsafe relative path in manifest: ' + row[field])
        archived = root / row['archive_path']
        if root not in archived.resolve().parents:
            raise SystemExit('Archive path leaves repository: ' + row['archive_path'])
        if row['archive_path'] not in verified:
            if (not archived.is_file() or archived.stat().st_size != int(row['bytes'])
                    or sha256(archived) != row['sha256']):
                raise SystemExit('Archive checksum mismatch: ' + str(archived))
            verified.add(row['archive_path'])
    dest.mkdir(parents=True, exist_ok=True)
    restored = set()
    adapted = []
    for row in rows:
        rel = Path(row['source_path'])
        if rel.is_absolute() or '..' in rel.parts:
            raise SystemExit('Unsafe relative path in manifest: ' + str(rel))
        if str(rel) in restored:
            continue
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / row['archive_path'], target)
        if row['role'] == 'code' and original_prefix:
            text = target.read_text(encoding='utf-8')
            changed = text.replace(original_prefix, str(dest))
            if text != changed:
                target.write_text(changed, encoding='utf-8')
                adapted.append(str(rel))
        restored.add(str(rel))
    # Verification before adaptation was done for every physical archive file.
    report = {'restored_files': len(restored), 'root_prefix_adapted_scripts': adapted,
              'archive': str(root), 'output': str(dest)}
    (dest / 'RESTORE_REPORT.json').write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))
    print('Run the commands in provenance/entrypoints.json with this directory as project root.')
    if not original_prefix:
        print('No absolute path prefix was supplied; copied code was not rewritten.')


if __name__ == '__main__':
    main()
