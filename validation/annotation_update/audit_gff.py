"""Read-only audit of explicit GFF3 Parent relationships; does not repair or infer them.

Usage: python3 audit_gff.py annotation.gff > audit.json
This checks annotation structure, not FASTA sequence identity or simulator truth.
"""
import collections
import hashlib
import json
import sys
from pathlib import Path
from urllib.parse import unquote


def records(path):
    with path.open() as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith('#'):
                continue
            columns = line.rstrip('\r\n').split('\t')
            if len(columns) != 9:
                raise ValueError(f'Line {number}: expected nine columns')
            attrs = dict(field.split('=', 1) for field in columns[8].split(';') if '=' in field)
            yield number, columns, attrs


def audit(path):
    needed = set()
    counts = collections.Counter()
    for _, columns, attrs in records(path):
        counts[columns[2]] += 1
        needed.update(unquote(v) for v in attrs.get('Parent', '').split(',') if v)
    parents = {}
    for number, columns, attrs in records(path):
        identifier = unquote(attrs.get('ID', ''))
        if identifier in needed:
            parents.setdefault(identifier, []).append({
                'line': number, 'id': identifier, 'reference': columns[0], 'type': columns[2],
                'start': int(columns[3]), 'end': int(columns[4]), 'strand': columns[6]})
    issues = []
    affected = set()
    for number, columns, attrs in records(path):
        for raw_parent in attrs.get('Parent', '').split(','):
            if not raw_parent:
                continue
            parent_id = unquote(raw_parent)
            matches = parents.get(parent_id, [])
            child = {'line': number, 'id': unquote(attrs.get('ID', '')), 'reference': columns[0],
                     'type': columns[2], 'start': int(columns[3]), 'end': int(columns[4]),
                     'strand': columns[6]}
            reason = None
            if not matches:
                reason = 'missing_parent'
            elif len(matches) != 1:
                # Legal discontinuous GFF3 parents need a richer model than this audit.
                reason = 'multiple_parent_rows_require_review'
            else:
                parent = matches[0]
                if parent['reference'] != child['reference']:
                    reason = 'different_references'
                elif parent['strand'] in '+-' and child['strand'] in '+-' and parent['strand'] != child['strand']:
                    reason = 'incompatible_strands'
                elif child['start'] < parent['start'] or child['end'] > parent['end']:
                    reason = 'outside_parent'
            if reason:
                affected.add(parent_id)
                issues.append({'reason': reason, 'child': child, 'parent_id': parent_id, 'parents': matches})
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'path': str(path.resolve()), 'sha256': digest.hexdigest(), 'bytes': path.stat().st_size,
            'records': sum(counts.values()), 'feature_types': dict(counts),
            'issue_counts': dict(collections.Counter(issue['reason'] for issue in issues)),
            'affected_parent_count': len(affected), 'issues': issues}


if __name__ == '__main__':
    json.dump(audit(Path(sys.argv[1])), sys.stdout, indent=2)
    print()
