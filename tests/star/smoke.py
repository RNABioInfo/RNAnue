"""Exercise the real pinned executable, independently of RNAnue's aligner."""
import argparse
import os
from pathlib import Path
import random
import subprocess
import tempfile


def run(args, **kwargs):
    result = subprocess.run([str(a) for a in args], capture_output=True,
                            text=True, timeout=90, **kwargs)
    if result.returncode:
        raise RuntimeError(f'Command failed ({result.returncode}): {args}\n'
                           f'{result.stdout}\n{result.stderr}')
    return result


def smoke(star, locator=None):
    star = Path(star).resolve()
    assert run([star, '--version']).stdout.strip() == '2.7.11b'
    with tempfile.TemporaryDirectory(prefix='rnanue STAR smoke ') as temporary:
        root = Path(temporary)
        if locator:
            resolved = Path(run([locator], cwd=root).stdout.strip()).resolve()
            assert resolved == star, (resolved, star)
        rng = random.Random(2411)
        sequence = ''.join(rng.choice('ACGT') for _ in range(4000))
        genome = root / 'reference.fa'
        genome.write_text('>chrTest\n' + sequence + '\n')
        reads = root / 'reads.fastq'
        reads.write_text(''.join(f'@read{i}\n{sequence[p:p+80]}\n+\n' + 'I'*80 + '\n'
                                 for i, p in enumerate([100, 800, 1600])))
        index = root / 'index'
        index.mkdir()
        run([star, '--runMode', 'genomeGenerate', '--genomeDir', index,
             '--genomeFastaFiles', genome, '--genomeSAindexNbases', '3',
             '--genomeChrBinNbits', '12', '--runThreadN', '2',
             '--outFileNamePrefix', str(root / 'index-log.')])
        run([star, '--genomeDir', index, '--readFilesIn', reads, '--runThreadN', '2',
             '--outSAMtype', 'SAM', '--outFileNamePrefix', str(root / 'mapping.')])
        lines = (root / 'mapping.Aligned.out.sam').read_text().splitlines()
        records = [line.split('\t') for line in lines if not line.startswith('@')]
        assert len(records) == 3, records
        assert {r[0] for r in records} == {'read0', 'read1', 'read2'}
        assert sorted(int(r[3]) for r in records) == [101, 801, 1601]
        assert all(r[2] == 'chrTest' and r[5] == '80M' and not int(r[1]) & 4 for r in records)
    print('STAR version, discovery, index generation and exact mapping passed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--star', required=True)
    parser.add_argument('--locator')
    args = parser.parse_args()
    smoke(args.star, args.locator)
