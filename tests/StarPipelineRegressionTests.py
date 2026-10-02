"""Real bundled STAR, cache lifecycle and unchanged downstream output contracts."""

import gzip
import hashlib
import os
import random
import shutil
import signal
import struct
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

BINARY = Path(sys.argv.pop(1)).resolve()
ENABLED = sys.argv.pop(1) == "ON"


def bam(path):
    """Small independent BAM reader; no optional Python packages or samtools needed."""
    data = gzip.open(path, "rb").read()
    assert data[:4] == b"BAM\1"
    size = struct.unpack_from("<i", data, 4)[0]
    header = data[8 : 8 + size].decode()
    offset = 8 + size
    refs = []
    count = struct.unpack_from("<i", data, offset)[0]
    offset += 4
    for _ in range(count):
        size = struct.unpack_from("<i", data, offset)[0]
        offset += 4
        refs.append(data[offset : offset + size - 1].decode())
        offset += size + 4
    records = []
    while offset < len(data):
        size = struct.unpack_from("<i", data, offset)[0]
        offset += 4
        r = data[offset : offset + size]
        offset += size
        ref, pos, packed, flags, length, mate_ref, mate_pos, tlen = struct.unpack_from(
            "<iiIIiiii", r
        )
        name_length, cigar_count = packed & 255, flags & 65535
        qname = r[32 : 32 + name_length - 1].decode()
        start = 32 + name_length
        cigar = "".join(
            str(c >> 4) + "MIDNSHP=X"[c & 15]
            for c in struct.unpack_from("<" + "I" * cigar_count, r, start)
        )
        start += cigar_count * 4
        bases = "=ACMGRSVTWYHKDBN"
        sequence = "".join(
            bases[(r[start + i // 2] >> (0 if i % 2 else 4)) & 15]
            for i in range(length)
        )
        start += (length + 1) // 2
        qualities = r[start : start + length]
        start += length
        tags = {}
        types = {
            "c": "b",
            "C": "B",
            "s": "h",
            "S": "H",
            "i": "i",
            "I": "I",
            "f": "f",
            "d": "d",
        }
        while start < len(r):
            name, kind = r[start : start + 2].decode(), chr(r[start + 2])
            start += 3
            if kind in "ZH":
                end = r.index(0, start)
                value = r[start:end].decode()
                start = end + 1
            elif kind == "A":
                value = chr(r[start])
                start += 1
            elif kind == "B":
                subtype = chr(r[start])
                n = struct.unpack_from("<i", r, start + 1)[0]
                start += 5
                value = struct.unpack_from("<" + types[subtype] * n, r, start)
                start += struct.calcsize(types[subtype]) * n
            else:
                value = struct.unpack_from("<" + types[kind], r, start)[0]
                start += struct.calcsize(types[kind])
            tags[name] = value
        records.append(
            {
                "name": qname,
                "ref": refs[ref] if ref >= 0 else None,
                "pos": pos,
                "flag": flags >> 16,
                "cigar": cigar,
                "sequence": sequence,
                "qualities": qualities,
                "mate_ref": refs[mate_ref] if mate_ref >= 0 else None,
                "mate_pos": mate_pos,
                "tlen": tlen,
                "tags": tags,
            }
        )
    return header, records


def fastq(name, sequence, quality="I"):
    return f"@{name}\n{sequence}\n+\n{quality * len(sequence)}\n"


class StarPipelineRegressions(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="RNAnue STAR pipeline ")
        self.root = Path(self.temp.name)
        rng = random.Random(2481)
        self.a, self.b = [
            "".join(rng.choice("ACGT") for _ in range(4000)) for _ in range(2)
        ]
        self.reference = self.root / "reference.fa"
        self.reference.write_text(f">chr1\n{self.a}\n>chr2\n{self.b}\n")
        self.index = self.root / "reference.fa.star_index"
        self.annotation = self.root / "features.gff"
        self.annotation.write_text(
            "##gff-version 3\n"
            + "".join(
                f"chr{i}\ttest\ttranscript\t1\t4000\t.\t+\t.\tID=T{i}\n" for i in [1, 2]
            )
        )
        self.input = self.root / "input" / "sample"
        self.input.mkdir(parents=True)
        self.output = self.root / "output with spaces"
        self.records = (
            fastq("chim", self.a[100:140] + self.b[800:840])
            + fastq("bg1", self.a[1200:1280])
            + fastq("bg2", self.b[2200:2280])
        )
        (self.input / "sample.fastq").write_text(self.records)

    def tearDown(self):
        self.temp.cleanup()

    def args(self, subcall, extra=()):
        return [
            subcall,
            "-T",
            str(self.input.parent),
            "-o",
            str(self.output),
            "-f",
            str(self.annotation),
            "-r",
            str(self.reference),
            "-t",
            "2",
            "--mask_multicopy_genes=false",
            "--deduplicate=false",
            "--trim_poly_g=false",
            "-q",
            "0",
            "-l",
            "1",
            "--min_complementarity=0",
            "--min_site_length_ratio=0",
            "--min_detect_length=1",
            "--orientation=both",
        ] + list(extra)

    def run_cli(self, subcall, extra=(), binary=BINARY, success=True):
        result = subprocess.run(
            [str(binary)] + self.args(subcall, extra),
            capture_output=True,
            text=True,
            timeout=90,
        )
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout + result.stderr

    def prepared(self, text=None):
        path = self.output / "01_preprocess/treatment/sample"
        path.mkdir(parents=True, exist_ok=True)
        with gzip.open(path / "sample_passed.fastq.gz", "wt") as f:
            f.write(self.records if text is None else text)
        return path

    def aligned(self):
        return self.output / "02_align/treatment/sample/sample_complete_aligned.bam"

    def assert_clean(self):
        bad = [
            p
            for p in self.root.rglob("*")
            if p.name.startswith(".rnanue-star-")
            or p.name.endswith(".lock")
            or p.name.startswith("Log.")
            or p.name.endswith("SJ.out.tab")
        ]
        self.assertEqual(bad, [])

    def test_availability_before_any_pipeline_output(self):
        if not ENABLED:
            for subcall in ["align", "complete"]:
                self.assertIn(
                    "RNANUE_BUILD_STAR=OFF",
                    self.run_cli(subcall, ["-a", "star"], success=False),
                )
                self.assertFalse(self.output.exists())
            return
        binary = self.root / "install/bin/RNAnue"
        binary.parent.mkdir(parents=True)
        shutil.copy2(BINARY, binary)
        for subcall in ["align", "complete"]:
            self.assertIn(
                "missing or not executable",
                self.run_cli(subcall, ["-a", "star"], binary, False),
            )
            self.assertFalse(self.output.exists())
        star = self.root / "install/libexec/rnanue/STAR"
        star.parent.mkdir(parents=True)
        star.write_text('#!/bin/sh\nprintf "2.7.10a\\n"\n')
        star.chmod(0o755)
        self.assertIn(
            "version mismatch", self.run_cli("complete", ["-a", "star"], binary, False)
        )
        self.assertFalse(self.output.exists())
        self.assert_clean()

    def test_segemehl_execution_and_legacy_empty_detection(self):
        self.prepared(fastq("ordinary", self.a[1200:1280]))
        self.run_cli("align")
        header, records = bam(self.aligned())
        self.assertIn("backend=segemehl", header)
        self.assertTrue(records)
        self.assertEqual({r["name"] for r in records}, {"ordinary"})
        self.run_cli("detect")
        self.assertFalse(self.index.exists())
        self.prepared("")
        self.run_cli("align")
        self.assertEqual(bam(self.aligned())[1], [])
        self.run_cli("detect")
        self.assert_clean()

    @unittest.skipUnless(ENABLED, "STAR-disabled build")
    def test_complete_and_standalone_detection_keep_output_contract(self):
        self.run_cli("complete", ["-a", "star"])
        header, records = bam(self.aligned())
        self.assertIn("RNAnue alignment backend=star;adapter=1;version=2.7.11b", header)
        self.assertEqual(len(records), 4)
        self.assertEqual({r["name"] for r in records}, {"chim", "bg1", "bg2"})
        splits = self.output / "03_detect/treatment/sample/sample_splits.bam"
        _, interactions = bam(splits)
        self.assertEqual(len(interactions), 2)
        self.assertEqual([r["name"] for r in interactions], ["chim", "chim"])
        self.assertEqual(
            {(r["ref"], r["pos"]) for r in interactions}, {("chr1", 100), ("chr2", 800)}
        )
        self.assertEqual(
            {r["sequence"] for r in interactions}, {self.a[100:140], self.b[800:840]}
        )
        for r in interactions:
            self.assertTrue({"XC", "XE", "XO", "XB"} <= r["tags"].keys())
            self.assertEqual(r["tags"]["XB"], 1)
        counts = (splits.parent / "sample_contiguous_transcript_counts.tsv").read_text()
        self.assertEqual(
            dict(line.split("\t") for line in counts.splitlines()),
            {"T1": "1", "T2": "1"},
        )
        self.assertTrue((self.output / "04_analyze").is_dir())
        self.assertTrue((self.output / "05_postprocess").is_dir())
        before = sorted(
            str(p.relative_to(self.output))
            for p in self.output.rglob("*")
            if p.is_file()
        )
        # No aligner option is required when detecting an existing RNAnue STAR BAM.
        self.run_cli("detect")
        self.assertEqual(bam(splits)[1], interactions)
        self.assertEqual(
            before,
            sorted(
                str(p.relative_to(self.output))
                for p in self.output.rglob("*")
                if p.is_file()
            ),
        )
        self.assertEqual(
            [p.name for p in self.aligned().parent.iterdir()], [self.aligned().name]
        )
        self.assert_clean()

    @unittest.skipUnless(ENABLED, "STAR-disabled build")
    def test_index_reuse_invalidation_lock_and_failed_rebuild(self):
        self.prepared()
        self.run_cli("align", ["-a", "star"])
        sa = self.index / "SA"
        before = sa.stat().st_mtime_ns
        self.run_cli("align", ["-a", "star", "--star_max_segment_gap=4"])
        self.assertEqual(sa.stat().st_mtime_ns, before)
        meta = self.index / "rnanue-index.meta"
        saved = meta.read_text()
        meta.write_text("incomplete")
        lock = self.index.with_name(self.index.name + ".lock")
        lock.mkdir()
        self.assertIn("locked", self.run_cli("align", ["-a", "star"], success=False))
        self.assertEqual(sa.stat().st_mtime_ns, before)
        lock.rmdir()
        self.run_cli("align", ["-a", "star"])
        self.assertEqual(meta.read_text(), saved)
        self.assertNotEqual(sa.stat().st_mtime_ns, before)
        (self.index / "SAindex").unlink()
        self.run_cli("align", ["-a", "star"])
        self.assertGreater((self.index / "SAindex").stat().st_size, 0)
        self.reference.write_text(
            self.reference.read_text().replace(self.a, "C" + self.a[1:])
        )
        self.run_cli("align", ["-a", "star"])
        self.assertNotEqual(meta.read_text(), saved)
        self.assert_clean()

    @unittest.skipUnless(ENABLED, "STAR-disabled build")
    def test_empty_short_and_malformed_compressed_inputs(self):
        folder = self.prepared("")
        self.run_cli("align", ["-a", "star"])
        header, records = bam(self.aligned())
        self.assertIn("backend=star", header)
        self.assertEqual(records, [])
        self.run_cli("detect")  # legacy empty and STAR empty inputs are both supported.
        with gzip.open(folder / "sample_passed.fastq.gz", "wt") as f:
            f.write(
                fastq("short", self.a[1200:1229]) + fastq("eligible", self.a[1200:1230])
            )
        self.run_cli("align", ["-a", "star", "--min_alignment_length=30"])
        self.assertEqual({r["name"] for r in bam(self.aligned())[1]}, {"eligible"})
        previous = self.aligned().read_bytes()
        with gzip.open(folder / "sample_passed.fastq.gz", "wt") as f:
            f.write("@broken\nAAA\n+\nI\n")
        self.assertIn(
            "Malformed FASTQ", self.run_cli("align", ["-a", "star"], success=False)
        )
        self.assertEqual(self.aligned().read_bytes(), previous)
        self.assert_clean()

    @unittest.skipUnless(ENABLED, "STAR-disabled build")
    def test_paired_categories_and_short_mate_remain_unsupported(self):
        folder = self.prepared()
        (folder / "sample_passed.fastq.gz").unlink()
        content = {
            "_merged_passed.fastq.gz": fastq(
                "merged", self.a[100:140] + self.b[800:840]
            ),
            "_singleton_passed_R1.fastq.gz": fastq("singleF", self.a[1200:1280]),
            "_singleton_passed_R2.fastq.gz": "",
            "_non_merged_passed_R1.fastq.gz": fastq("pair/1", self.a[1200:1280]),
            "_non_merged_passed_R2.fastq.gz": fastq("pair/2", self.b[2200:2220]),
        }
        for suffix, text in content.items():
            with gzip.open(folder / ("sample" + suffix), "wt") as f:
                f.write(text)
        self.run_cli("align", ["-a", "star", "--min_alignment_length=30"])
        header, records = bam(self.aligned())
        self.assertIn("backend=star", header)
        pairs = [r for r in records if r["name"].startswith("pair")]
        self.assertTrue(pairs)
        self.assertTrue(all(r["flag"] & 1 for r in pairs))
        self.assertTrue(all(len(r["sequence"]) >= 30 for r in pairs))
        self.run_cli("detect")
        _, interactions = bam(
            self.output / "03_detect/treatment/sample/sample_splits.bam"
        )
        self.assertEqual([r["name"] for r in interactions], ["merged", "merged"])
        # Unequal mate counts fail before launching STAR.
        with gzip.open(folder / "sample_non_merged_passed_R2.fastq.gz", "wt") as f:
            f.write("")
        self.assertIn(
            "unequal record counts",
            self.run_cli("align", ["-a", "star"], success=False),
        )
        self.assert_clean()

    @unittest.skipUnless(ENABLED, "STAR-disabled build")
    def test_coverage_rounding_and_effective_argument_profile(self):
        self.prepared(fastq("clipped", self.a[1200:1280] + "N" * 21))
        self.run_cli(
            "align",
            [
                "-a",
                "star",
                "--min_split_coverage=80",
                "--allow_multimapping=false",
                "--star_max_multimaps=17",
                "--min_fragment_length=16",
            ],
        )
        header, records = bam(self.aligned())
        self.assertEqual(len(records), 1)
        self.assertIn("80M21S", records[0]["cigar"])
        for name, value in [
            ("outFilterMultimapNmax", 1),
            ("chimMultimapNmax", 1),
            ("chimJunctionOverhangMin", 16),
        ]:
            self.assertIn(f'"--{name}" "{value}"', header)
        self.run_cli("align", ["-a", "star", "--min_split_coverage=81"])
        self.assertEqual(bam(self.aligned())[1], [])
        self.assert_clean()

    @unittest.skipUnless(ENABLED, "STAR-disabled build")
    def test_ordinary_intron_alignment_is_decoded_as_two_segments(self):
        self.a = self.a[:140] + "GT" + self.a[142:188] + "AG" + self.a[190:]
        self.reference.write_text(f">chr1\n{self.a}\n>chr2\n{self.b}\n")
        self.prepared(fastq("spliced", self.a[100:140] + self.a[190:230]))
        self.run_cli("align", ["-a", "star", "--star_max_intron_length=100"])
        _, records = bam(self.aligned())
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["cigar"], "40M50N40M")
        self.run_cli("detect")
        _, interactions = bam(
            self.output / "03_detect/treatment/sample/sample_splits.bam"
        )
        self.assertEqual(len(interactions), 2)
        self.assertEqual({r["pos"] for r in interactions}, {100, 190})
        self.assertEqual(
            {r["sequence"] for r in interactions}, {self.a[100:140], self.a[190:230]}
        )
        self.assert_clean()

    @unittest.skipUnless(ENABLED, "STAR-disabled build")
    def test_failed_index_process_and_interruption_preserve_previous_index(self):
        self.prepared()
        self.run_cli("align", ["-a", "star"])
        saved = (self.index / "SA").read_bytes()
        saved_metadata = (self.index / "rnanue-index.meta").read_bytes()
        self.reference.write_text(
            self.reference.read_text().replace(self.a, "C" + self.a[1:])
        )
        binary = self.root / "install/bin/RNAnue"
        binary.parent.mkdir(parents=True)
        shutil.copy2(BINARY, binary)
        star = self.root / "install/libexec/rnanue/STAR"
        star.parent.mkdir(parents=True)
        preamble = '#!/bin/sh\nif [ "$1" = "--version" ]; then printf "2.7.11b\\n"; exit 0; fi\n'
        star.write_text(preamble + 'echo "intentional STAR failure" >&2\nexit 17\n')
        star.chmod(0o755)
        text = self.run_cli("align", ["-a", "star"], binary, False)
        self.assertIn("exit status 17", text)
        self.assertIn("intentional STAR failure", text)
        self.assertEqual((self.index / "SA").read_bytes(), saved)
        self.assertEqual(
            (self.index / "rnanue-index.meta").read_bytes(), saved_metadata
        )
        self.assert_clean()
        star.write_text(preamble + "exec /bin/sleep 30\n")
        process = subprocess.Popen(
            [str(binary)] + self.args("align", ["-a", "star"]),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        try:
            deadline = time.monotonic() + 10
            while (
                not list(self.root.glob(".rnanue-star-*/process.log"))
                and time.monotonic() < deadline
            ):
                if process.poll() is not None:
                    break
                time.sleep(0.02)
            self.assertTrue(list(self.root.glob(".rnanue-star-*/process.log")))
            process.send_signal(signal.SIGTERM)
            stdout, stderr = process.communicate(timeout=10)
            self.assertNotEqual(process.returncode, 0)
            self.assertIn("interrupted by signal", stdout + stderr)
        finally:
            if process.poll() is None:
                process.kill()
                process.communicate()
        self.assertEqual((self.index / "SA").read_bytes(), saved)
        self.assert_clean()

    @unittest.skipUnless(
        ENABLED and os.geteuid() != 0, "requires STAR and non-root permissions"
    )
    def test_read_only_reference_directory_uses_output_index(self):
        self.prepared()
        directory = self.root / "read only reference"
        directory.mkdir()
        self.reference.rename(directory / self.reference.name)
        self.reference = directory / self.reference.name
        directory.chmod(0o555)
        try:
            self.run_cli("align", ["-a", "star"])
            index = self.output / "reference.fa.star_index"
            self.assertTrue((index / "rnanue-index.meta").is_file())
            self.assertFalse((directory / "reference.fa.star_index").exists())
            before = (index / "SA").stat().st_mtime_ns
            self.run_cli("align", ["-a", "star"])
            self.assertEqual((index / "SA").stat().st_mtime_ns, before)
        finally:
            directory.chmod(0o755)
        self.assert_clean()

    def exercise_supplied_index(self, backend):
        ordinary = fastq("ordinary", self.a[1200:1280])
        self.prepared(ordinary)
        (self.input / "sample.fastq").write_text(ordinary)
        args = ["-a", backend]
        self.run_cli("align", args)
        expected = bam(self.aligned())[1]
        expected_files = sorted(self.aligned().parent.iterdir())
        automatic = self.index if backend == "star" else self.reference.with_suffix(".idx")
        store = self.root / "precomputed elsewhere"
        store.mkdir()
        supplied = store / "arbitrary index name"
        automatic.rename(supplied)
        if backend == "star":
            (supplied / "rnanue-index.meta").unlink()
            automatic.mkdir()
            poison = automatic / "do not replace"
        else:
            poison = automatic
        poison.write_text("Explicit index must take precedence")
        link = self.root / "index symlink"
        link.symlink_to(supplied, target_is_directory=backend == "star")

        def snapshot():
            paths = [supplied] + (sorted(supplied.rglob("*")) if supplied.is_dir() else [])
            return [(str(p.relative_to(store)), p.stat().st_mtime_ns, p.stat().st_mode,
                     hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None)
                    for p in paths]

        def readonly(enabled):
            paths = [supplied] + (list(supplied.rglob("*")) if supplied.is_dir() else [])
            for path in paths:
                path.chmod((0o555 if enabled else 0o755) if path.is_dir()
                           else (0o444 if enabled else 0o644))

        for form in ["absolute", "relative symlink", "config"]:
            if backend == "star" and form != "absolute":
                (supplied / "rnanue-index.meta").write_text("unrelated RNAnue metadata")
                # A different reference fingerprint must not invalidate an explicit index.
                first = "C" if self.a[0] != "C" else "A"
                self.reference.write_text(f">chr1\n{first}{self.a[1:]}\n>chr2\n{self.b}\n")
            readonly(True)
            before = snapshot()
            try:
                if form == "absolute":
                    options = ["--alignment_index", str(supplied)]
                elif form == "relative symlink":
                    options = ["-i", os.path.relpath(link)]
                else:
                    config = self.root / "index configuration.cfg"
                    config.write_text(f"alignment_index = {link}\nmask_multicopy_genes = false\n")
                    options = ["-c", str(config)]
                text = self.run_cli("align", args + options)
                self.assertIn("Using supplied", text)
                self.assertNotIn("Building", text)
                self.assertEqual(bam(self.aligned())[1], expected)
                self.assertEqual(snapshot(), before)
                self.assertEqual(poison.read_text(), "Explicit index must take precedence")
                self.assertEqual(list(store.iterdir()), [supplied])
                self.assertEqual(sorted(self.aligned().parent.iterdir()), expected_files)
                self.assert_clean()
            finally:
                readonly(False)

        readonly(True)
        before = snapshot()
        try:
            self.run_cli("complete", args + ["-i", str(supplied)])
            self.assertEqual(bam(self.aligned())[1], expected)
            self.assertEqual(snapshot(), before)
            self.assertEqual(poison.read_text(), "Explicit index must take precedence")
            self.assert_clean()
        finally:
            readonly(False)

    def test_supplied_segemehl_index(self):
        self.exercise_supplied_index("segemehl")

    @unittest.skipUnless(ENABLED, "STAR-disabled build")
    def test_supplied_star_index(self):
        self.exercise_supplied_index("star")

    @unittest.skipUnless(ENABLED, "STAR-disabled build")
    def test_malformed_supplied_star_index_does_not_fall_back(self):
        self.prepared(fastq("ordinary", self.a[1200:1280]))
        self.run_cli("align", ["-a", "star"])
        existing = (self.index / "SA").stat().st_mtime_ns
        previous = self.aligned().read_bytes()
        malformed = self.root / "empty provided index"
        malformed.mkdir()
        text = self.run_cli("align", ["-a", "star", "-i", str(malformed)], success=False)
        self.assertIn("STAR failed", text)
        self.assertEqual(list(malformed.iterdir()), [])
        self.assertEqual((self.index / "SA").stat().st_mtime_ns, existing)
        self.assertEqual(self.aligned().read_bytes(), previous)
        self.assert_clean()


if __name__ == "__main__":
    unittest.main()
