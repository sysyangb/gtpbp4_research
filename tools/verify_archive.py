#!/usr/bin/env python3
"""Verify retained files, source mappings, and generated repository metadata."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def relative_file(root, value):
    rel = Path(value)
    if not value or rel.is_absolute() or '..' in rel.parts:
        raise ValueError('Unsafe relative file path: ' + value)
    target = root / rel
    if root.resolve() not in target.resolve().parents:
        raise ValueError('File path leaves repository: ' + value)
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path,
                        help='Also verify original files against their source checksums')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    with (root / 'provenance/file_manifest.tsv').open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    failures, checked, source_checked = [], {}, {}
    for row in rows:
        try:
            path = relative_file(root, row['archive_path'])
        except ValueError as error:
            failures.append(str(error))
            continue
        if str(path) not in checked:
            if not path.is_file():
                failures.append('Missing: ' + str(path))
                continue
            checked[str(path)] = digest(path)
        if checked[str(path)] != row['sha256'] or path.stat().st_size != int(row['bytes']):
            failures.append('Archive changed: ' + row['archive_path'])
        if args.source_root:
            try:
                source = relative_file(args.source_root, row['source_path'])
            except ValueError as error:
                failures.append(str(error))
                continue
            expected_hash = row.get('source_sha256') or row['sha256']
            expected_bytes = int(row.get('source_bytes') or row['bytes'])
            if source.is_file() and str(source) not in source_checked:
                source_checked[str(source)] = digest(source)
            if (not source.is_file() or source.stat().st_size != expected_bytes
                    or source_checked.get(str(source)) != expected_hash):
                failures.append('Original changed or missing: ' + row['source_path'])
    catalog = json.loads((root / 'provenance/panel_catalog.json').read_text())
    mappings = {(row['source_path'], row['archive_path']): row for row in rows}
    for panel in catalog:
        for item in panel['files']:
            mapped = mappings.get((item['source_path'], item['archive_path']))
            if (not mapped or item['sha256'] != mapped['sha256']
                    or int(item['bytes']) != int(mapped['bytes'])):
                failures.append('Broken panel mapping: ' + panel['panel'] + '/' + item['archive_path'])
    generated_file = root / 'provenance/generated_manifest.tsv'
    generated_count = 0
    if generated_file.exists():
        with generated_file.open(newline='') as f:
            generated = list(csv.DictReader(f, delimiter='\t'))
        for item in generated:
            try:
                saved = relative_file(root, item['archive_path'])
            except ValueError as error:
                failures.append(str(error))
                continue
            if (not saved.is_file() or digest(saved) != item['sha256']
                    or ('bytes' in item and saved.stat().st_size != int(item['bytes']))):
                failures.append('Generated documentation or tool changed: ' + item['archive_path'])
        generated_count = len(generated)
    report = {'ok': not failures, 'source_mappings': len(rows),
              'physical_source_files': len(checked), 'panels': len(catalog),
              'physical_source_bytes': sum(Path(p).stat().st_size for p in checked),
              'generated_files_verified': generated_count,
              'originals_verified': bool(args.source_root), 'failures': failures}
    if args.report:
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))
    raise SystemExit(0 if report['ok'] else 1)


if __name__ == '__main__':
    main()
