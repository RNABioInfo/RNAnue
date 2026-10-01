"""CLI integrity regressions; each case uses a fresh directory, never resumption."""
import gzip
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

BINARY = Path(sys.argv.pop(1)).resolve()


def fastq(name, sequence='ACGT' * 10, quality=None):
    return f'@{name}\n{sequence}\n+\n{quality if quality is not None else "I" * len(sequence)}\n'


class PipelineCliRegressions(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='rnanue-cli-')
        self.root = Path(self.temp.name)
        self.input = self.root/'input'/'sample'
        self.input.mkdir(parents=True)
        self.annotation = self.root/'features.gff'
        self.annotation.write_text('##gff-version 3\nchr1\ttest\ttranscript\t1\t1000\t.\t+\t.\tID=T\n')

    def tearDown(self):
        self.temp.cleanup()

    def run_case(self, forward, reverse=None, deduplicate=False, threads=2, extra=(), output_name='out',
                 min_read_length=1):
        (self.input/('sample_R1.fastq' if reverse is not None else 'sample.fastq')).write_text(forward)
        if reverse is not None:
            (self.input/'sample_R2.fastq').write_text(reverse)
        command = [str(BINARY),'preprocess','-T',str(self.input.parent),'-o',str(self.root/output_name),
                   '-f',str(self.annotation),'--threads',str(threads),'--min_read_quality','0',
                   '--min_read_length',str(min_read_length),
                   '--no_trim_poly_g']
        if not deduplicate:
            command.append('--no_deduplicate')
        return subprocess.run(command+list(extra),capture_output=True,text=True,timeout=30)

    def output_count(self):
        return sum(len(gzip.open(p,'rt').read().splitlines())//4 for p in (self.root/'out').rglob('*.fastq.gz'))

    def output_records(self, output_name='out'):
        records = []
        for path in (self.root/output_name).rglob('*.fastq.gz'):
            with gzip.open(path, 'rt') as handle:
                lines = handle.read().splitlines()
            records.extend((lines[i], lines[i+1], lines[i+3]) for i in range(0, len(lines), 4))
        return sorted(records)

    def test_unequal_mates_fail_clearly(self):
        result = self.run_case(fastq('one/1')+fastq('two/1'),fastq('one/2'))
        self.assertEqual(result.returncode,1,result.stdout+result.stderr)
        self.assertIn('pair 2',result.stdout+result.stderr)
        self.assertIn('<EOF>',result.stdout+result.stderr)

    def test_unrelated_names_fail(self):
        result = self.run_case(fastq('one'),fastq('two'))
        self.assertEqual(result.returncode,1,result.stdout+result.stderr)
        self.assertIn('IDs: one / two',result.stdout+result.stderr)

    def test_common_names_succeed(self):
        result = self.run_case(fastq('one/1 1:N:0:ACGT'),fastq('one/2 2:N:0:ACGT'))
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual(self.output_count(),1)

    def test_invalid_pair_buffer_is_rejected(self):
        result = self.run_case(fastq('one'),fastq('one'),extra=['--chunk_size','1'])
        self.assertEqual(result.returncode,1,result.stdout+result.stderr)
        self.assertIn('--chunk_size >= 2',result.stdout+result.stderr)

    def test_malformed_fastq_is_controlled_error(self):
        result = self.run_case(fastq('one',quality='I'))
        self.assertEqual(result.returncode,1,result.stdout+result.stderr)
        self.assertIn('FASTQ record 1',result.stdout+result.stderr)

    def test_duplicate_id_does_not_retain_unselected_records(self):
        result = self.run_case(fastq('same')*3+fastq('same','AAAA'*10),deduplicate=True,threads=4)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual(self.output_count(),2)

    def test_repeated_adapters_single_end_and_dimers(self):
        adapter5 = 'TGTAGATCTCGGTGGTCGCCGTATCATT'
        adapter3 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
        insert = 'TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCTCTCACCG'
        sequence = adapter5*2 + insert + adapter3*2
        qualities = ''.join(chr(40+i % 30) for i in range(len(sequence)))
        result = self.run_case(fastq('keep', sequence, qualities) +
                               fastq('dimer', adapter5*2+adapter3*2),
                               min_read_length=0, extra=['--adapter_5p_forward', adapter5,
                                                        '--adapter_3p_forward', adapter3])
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertEqual(self.output_records(),
                         [('@keep', insert, qualities[len(adapter5)*2:len(adapter5)*2+len(insert)])])

    def test_repeated_adapters_on_both_paired_mates(self):
        adapter = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
        forward_insert, reverse_insert = 'C'*40, 'AC'*20
        result = self.run_case(fastq('keep/1', adapter*2+forward_insert) + fastq('dimer/1', adapter*2),
                               fastq('keep/2', reverse_insert+adapter*2) + fastq('dimer/2', adapter*2),
                               min_read_length=0, extra=['--adapter_5p_forward', adapter,
                                                        '--adapter_3p_reverse', adapter])
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertEqual(sorted(r[1] for r in self.output_records()), [reverse_insert, forward_insert])
        self.assertTrue(all(len(r[1]) == len(r[2]) for r in self.output_records()))

    def test_simulation_adapter_probes_and_thread_equivalence(self):
        source = ('TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCTCTCACCGCCGCGGCCCGGGTTCGATT'
                  'CCCGGTCAGGGAAAGA')
        adapter = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
        input_text, expected = '', []
        for length in (40, 60, 73, 75, 84, 90, 100, 120):
            insert = (source*2)[:length]
            sequence = (insert+adapter*5)[:150]
            input_text += fastq(f'insert_{length}', sequence)
            expected.append((f'@insert_{length}', insert, 'I'*length))
        for threads in (2, 4):
            name = f'threads_{threads}'
            result = self.run_case(input_text, threads=threads, output_name=name,
                                   extra=['--adapter_3p_forward', adapter])
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertEqual(self.output_records(name), sorted(expected))


if __name__ == '__main__':
    unittest.main()
