#!/usr/bin/env python3
"""Check distributed paths and Office packages for excluded document/media assets."""
import argparse
import gzip
import json
from pathlib import Path
import re
import subprocess
import zipfile

BLOCKED_EXTENSIONS = {
    '.png', '.jpg', '.jpeg', '.tif', '.tiff', '.svg', '.pdf', '.eps', '.ps',
    '.ai', '.emf', '.wmf', '.gif', '.webp', '.bmp', '.ico', '.heic', '.avif',
    '.doc', '.docx', '.ppt', '.pptx', '.rtf', '.odt', '.rmd', '.html', '.htm',
}
IMAGE_PREFIXES = (b'\x89PNG\r\n\x1a\n', b'\xff\xd8\xff', b'%PDF-',
                  b'II*\x00', b'MM\x00*', b'GIF87a', b'GIF89a')


def check(root, paths):
    failures, workbooks = [], []
    for rel in paths:
        p = root / rel
        if not p.is_file():
            failures.append({'path': rel, 'reason': 'missing file'})
            continue
        path = Path(rel)
        if path.name.startswith('._') or path.name == '.DS_Store':
            failures.append({'path': rel, 'reason': 'filesystem metadata'})
        if path.suffix.lower() in BLOCKED_EXTENSIONS or 'manuscript' in path.parts:
            failures.append({'path': rel, 'reason': 'document or image path'})
        if 'notes' in path.parts or re.search(r'(?:legend|caption|figure_contract|methods_note)', path.name, re.I):
            failures.append({'path': rel, 'reason': 'caption or narrative asset path'})
        with p.open('rb') as f:
            prefix = f.read(128)
        if prefix.startswith(IMAGE_PREFIXES) or prefix.lstrip().startswith(b'<svg'):
            failures.append({'path': rel, 'reason': 'image or PDF signature'})
        if p.suffix.lower() == '.gz':
            try:
                with gzip.open(p, 'rb') as f:
                    decoded = f.read(128)
                if decoded.startswith(IMAGE_PREFIXES) or decoded.lstrip().startswith(b'<svg'):
                    failures.append({'path': rel, 'reason': 'compressed image or PDF'})
            except (OSError, EOFError) as exc:
                failures.append({'path': rel, 'reason': 'invalid gzip: ' + str(exc)})
        if p.suffix.lower() == '.xlsx':
            try:
                with zipfile.ZipFile(p) as z:
                    media = [n for n in z.namelist() if re.search(r'(?:^|/)(?:media|drawings|charts|embeddings)/', n)]
                    anchors = []
                    for n in z.namelist():
                        if re.fullmatch(r'xl/worksheets/sheet\d+\.xml', n):
                            if re.search(rb'<(?:\w+:)?(?:drawing|legacyDrawing|picture|oleObjects)\b', z.read(n)):
                                anchors.append(n)
                    workbooks.append({'path': rel, 'media_parts': len(media), 'drawing_sheets': len(anchors)})
                    if media or anchors:
                        failures.append({'path': rel, 'reason': 'embedded graphic or object', 'parts': media + anchors})
            except (OSError, zipfile.BadZipFile) as exc:
                failures.append({'path': rel, 'reason': 'invalid workbook: ' + str(exc)})
    return {'ok': not failures, 'files_checked': len(paths), 'workbooks': workbooks, 'failures': failures}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tracked', action='store_true', help='Check only Git tracked paths')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    if args.tracked:
        raw = subprocess.check_output(['git', 'ls-files', '-z'], cwd=root)
        paths = [p.decode() for p in raw.split(b'\x00') if p]
    else:
        paths = sorted(str(p.relative_to(root)) for p in root.rglob('*')
                       if p.is_file() and '.git' not in p.relative_to(root).parts
                       and not p.name.startswith('._') and p.name != '.DS_Store')
    report = check(root, paths)
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report['ok'] else 1)


if __name__ == '__main__':
    main()
