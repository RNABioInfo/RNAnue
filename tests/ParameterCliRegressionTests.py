"""Public option spelling contract and side-effect-free backend rejection."""

import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BINARY = Path(sys.argv.pop(1)).resolve()
# Retired spellings are intentional: they must never be accepted as aliases or prefixes.
MIGRATION = {
    "trtms": "treatment_dir",
    "ctrls": "control_dir",
    "outdir": "output_dir",
    "loglevel": "log_level",
    "featuretypes": "feature_types",
    "maskmulticopy": "mask_multicopy_genes",
    "mincopyident": "min_copy_identity",
    "chunksize": "chunk_size",
    "trimpolyg": "trim_poly_g",
    "minpolygcount": "min_poly_g_count",
    "mtrim": "max_adapter_mismatch_rate",
    "minovltrim": "min_adapter_overlap",
    "minqual": "min_read_quality",
    "minlen": "min_read_length",
    "wqual": "min_window_quality",
    "wtrim": "quality_window_size",
    "adpt5f": "adapter_5p_forward",
    "adpt5r": "adapter_5p_reverse",
    "adpt3f": "adapter_3p_forward",
    "adpt3r": "adapter_3p_reverse",
    "minovl": "min_merge_overlap",
    "mmerge": "max_merge_mismatch_rate",
    "dbref": "reference_genome",
    "multimap": "allow_multimapping",
    "accuracy": "seg_alignment_accuracy",
    "minfragsco": "seg_min_fragment_score",
    "minalignlen": "min_alignment_length",
    "minfraglen": "min_fragment_length",
    "minsplicecov": "min_split_coverage",
    "maxprim": "max_primary_alignments",
    "mapqmin": "min_mapping_quality",
    "cmplmin": "min_complementarity",
    "sitelenratio": "min_site_length_ratio",
    "mindetectlen": "min_detect_length",
    "nrgmax": "max_hybridization_energy",
    "exclclipping": "exclude_soft_clipping",
    "splicing": "filter_splicing",
    "altsplice": "remove_alt_splicing",
    "spltol": "splicing_tolerance",
    "includewobble": "include_wobble",
    "mincontr": "min_hit_contribution",
    "maxselfoverlap": "max_self_overlap",
    "clustmethod": "clustering_strand_specificity",
    "clustdist": "clustering_distance",
    "clustfrac": "min_cluster_overlap",
    "padj": "max_adjusted_p_value",
    "mincount": "min_interaction_support",
    "mineffdens": "min_support_per_effective_bp",
    "maxcovcomp": "max_coverage_components",
    "minarmbal": "min_arm_balance",
    "intfrac": "min_interaction_overlap",
    "no-preprocess": "no_preprocess",
    "no-deduplicate": "no_deduplicate",
    "no-trimpolyg": "no_trim_poly_g",
    "no-maskmulticopy": "no_mask_multicopy_genes",
    "no-multimap": "no_multimapping",
    "keep-altsplice": "keep_alt_splicing",
}
STAR_CONTROLS = {
    "star_max_multimaps": 10,
    "star_min_junction_overhang": 15,
    "star_max_segment_gap": 3,
    "star_min_nonchimeric_score_drop": 10,
    "star_max_chimeric_score_drop": 30,
    "star_max_intron_length": 10,
}
BACKEND_MISMATCH = "requires aligner=segemehl"


class ParameterCliRegressions(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="rnanue-parameters-")
        self.root = Path(self.temp.name)
        self.input = self.root / "input"
        self.input.mkdir()
        # These would fail if complete reached annotation or preprocessing.
        self.annotation = self.root / "invalid.gff"
        self.annotation.write_text("deliberately malformed annotation\n")
        (self.input / "invalid.fastq").write_text("not fastq\n")
        self.output = self.root / "must not exist"
        self.config = self.root / "parameters with spaces.cfg"

    def tearDown(self):
        self.temp.cleanup()

    def run_cli(self, args):
        return subprocess.run(
            [str(BINARY)] + args,
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )

    def required(self, subcall="align"):
        return [
            subcall,
            "-T",
            str(self.input),
            "-o",
            str(self.output),
            "-f",
            str(self.annotation),
            "-r",
            str(self.root / "missing reference.fa"),
        ]

    def assert_failed(self, result, message):
        text = result.stdout + result.stderr
        self.assertNotEqual(result.returncode, 0, text)
        self.assertIn(message, text)
        self.assertFalse(self.output.exists(), text)

    def test_help_has_exact_canonical_names_and_shortcuts(self):
        result = self.run_cli(["--help"])
        self.assertEqual(result.returncode, 0, result.stderr)
        names = re.findall(
            r"^  (?:-[A-Za-z] \[ )?--([a-z][a-z0-9_]*)", result.stdout, re.MULTILINE
        )
        expected = (
            set(MIGRATION.values())
            | set(STAR_CONTROLS)
            | {
                "threads",
                "features",
                "orientation",
                "preprocess",
                "deduplicate",
                "config",
                "help",
                "version",
                "subcall",
                "aligner",
            }
        )
        self.assertEqual(set(names), expected)
        self.assertEqual(len(names), len(expected))
        shortcuts = dict(
            re.findall(
                r"^  -([A-Za-z]) \[ --([a-z][a-z0-9_]*) \]", result.stdout, re.MULTILINE
            )
        )
        self.assertEqual(
            shortcuts,
            {
                "T": "treatment_dir",
                "C": "control_dir",
                "t": "threads",
                "o": "output_dir",
                "r": "reference_genome",
                "f": "features",
                "a": "aligner",
                "c": "config",
                "q": "min_read_quality",
                "l": "min_read_length",
                "h": "help",
                "v": "version",
            },
        )
        self.assertIn("segemehl-specific", result.stdout)
        self.assertIn("STAR", result.stdout)

    def test_retired_long_names_and_abbreviations_fail(self):
        for name in list(MIGRATION) + [
            "alignment_accuracy",
            "min_fragment_score",
            "treat",
            "star_max_multi",
        ]:
            with self.subTest(option=name):
                self.assert_failed(
                    self.run_cli(["--help", "--" + name]), "unrecognised option"
                )

    def test_retired_config_keys_fail(self):
        for name in list(MIGRATION) + [
            "alignment_accuracy",
            "min_fragment_score",
            "treat",
        ]:
            with self.subTest(option=name):
                self.config.write_text(name + " = 1\n")
                self.assert_failed(
                    self.run_cli(self.required() + ["-c", str(self.config)]),
                    "unrecognised option",
                )

    def test_retired_shortcuts_fail(self):
        for flag in ["-p", "-s"]:
            self.assert_failed(
                self.run_cli(["--help", flag, "2"]), "unrecognised option"
            )
        self.assert_failed(
            self.run_cli(self.required() + ["-t", str(self.input)]), "threads"
        )
        for flag in ["-h", "-v"]:
            result = self.run_cli([flag])
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_seg_controls_with_star_rejected_before_pipeline_io_cli_and_config(self):
        for subcall in ["align", "complete"]:
            for config in [False, True]:
                with self.subTest(subcall=subcall, config=config):
                    args = self.required(subcall)
                    if config:
                        self.config.write_text("aligner = star\n")
                        args += ["-c", str(self.config)]
                    else:
                        args += ["-a", "star"]
                    self.assert_failed(
                        self.run_cli(args + ["--seg_alignment_accuracy=90"]),
                        BACKEND_MISMATCH,
                    )

    def test_explicit_star_controls_cannot_be_ignored_by_segemehl(self):
        for subcall in ["align", "complete"]:
            for name, value in STAR_CONTROLS.items():
                for config in [False, True]:
                    with self.subTest(subcall=subcall, option=name, config=config):
                        setting = f"{name}={value}"
                        args = self.required(subcall)
                        if config:
                            self.config.write_text(setting + "\n")
                            args += ["-c", str(self.config)]
                        else:
                            args += ["--" + setting]
                        self.assert_failed(self.run_cli(args), "requires aligner=star")

    def test_negative_zero_overflow_and_valid_star_boundaries(self):
        for name in STAR_CONTROLS:
            lower = (
                1 if name in ["star_max_multimaps", "star_min_junction_overhang"] else 0
            )
            for value in [-1, 2147483648, -2147483649] + ([0] if lower else []):
                for config in [False, True]:
                    with self.subTest(option=name, value=value, config=config):
                        args = self.required() + ["--aligner=star"]
                        if config:
                            self.config.write_text(f"{name} = {value}\n")
                            args += ["-c", str(self.config)]
                        else:
                            args += [f"--{name}={value}"]
                        result = self.run_cli(args)
                        self.assert_failed(result, name)
                        self.assertNotIn(
                            BACKEND_MISMATCH, result.stdout + result.stderr
                        )
            for value in [lower, 2147483647]:
                with self.subTest(option=name, valid=value):
                    self.assert_failed(
                        self.run_cli(
                            self.required()
                            + [
                                "--aligner=star",
                                "--seg_alignment_accuracy=90",
                                f"--{name}={value}",
                            ]
                        ),
                        BACKEND_MISMATCH,
                    )

    def test_invalid_backend_and_fragment_length(self):
        for backend in ["STAR", "Segemehl", "other"]:
            self.assert_failed(
                self.run_cli(self.required() + ["--aligner=" + backend]),
                "aligner must be exactly",
            )
        for length in ["0", "2147483648", "-1"]:
            self.assert_failed(
                self.run_cli(
                    self.required()
                    + ["--aligner=star", "--min_fragment_length=" + length]
                ),
                "min_fragment_length",
            )

    def test_cli_backend_precedence_is_preserved(self):
        self.config.write_text("aligner = segemehl\n")
        self.assert_failed(
            self.run_cli(
                self.required("complete")
                + [
                    "-c",
                    str(self.config),
                    "--aligner=star",
                    "--seg_alignment_accuracy=90",
                ]
            ),
            BACKEND_MISMATCH,
        )
        self.config.write_text("aligner = star\nstar_max_multimaps = 10\n")
        self.assert_failed(
            self.run_cli(
                self.required("complete")
                + ["-c", str(self.config), "--aligner=segemehl"]
            ),
            "requires aligner=star",
        )


if __name__ == "__main__":
    unittest.main()
