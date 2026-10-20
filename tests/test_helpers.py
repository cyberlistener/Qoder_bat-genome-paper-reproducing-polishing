#!/usr/bin/env python3
"""Minimal unit tests for the pure helper functions.

Run with:
    python -m unittest discover -s tests -v

These tests do not invoke HyPhy, R, g:Profiler or Snakemake.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from common import bh_fdr, codon_qc, parse_rate_distributions  # noqa: E402


class TestBH(unittest.TestCase):
    def test_known_values(self):
        q = bh_fdr([0.01, 0.02, 0.03, 0.04, 0.05])
        self.assertTrue(all(abs(x - 0.05) < 1e-12 for x in q))

    def test_nan_preserved_and_monotone(self):
        q = bh_fdr([0.001, float("nan"), 0.02, 0.5])
        self.assertTrue(math.isnan(q[1]))
        finite = q[[0, 2, 3]]
        self.assertTrue(all(finite[0] <= finite[1] <= finite[2]))
        self.assertTrue(all(0 <= x <= 1 for x in finite))


class TestRateDistributions(unittest.TestCase):
    def test_list_shape(self):
        omega, prop, n = parse_rate_distributions([[0.2, 0.9], [2.5, 0.1]])
        self.assertEqual(n, 2)
        self.assertAlmostEqual(omega, 2.5)
        self.assertAlmostEqual(prop, 0.1)

    def test_dict_shape(self):
        omega, prop, n = parse_rate_distributions({"0": [0.2, 0.9], "1": [3.0, 0.05]})
        self.assertEqual(n, 2)
        self.assertAlmostEqual(omega, 3.0)
        self.assertAlmostEqual(prop, 0.05)

    def test_missing(self):
        omega, prop, n = parse_rate_distributions(None)
        self.assertTrue(math.isnan(omega))
        self.assertEqual(prop, 0.0)
        self.assertEqual(n, 0)


class TestCodonQC(unittest.TestCase):
    def test_internal_stop(self):
        seqs = {"a": "TAAGCT", "b": "TAAGCT"}
        qc = codon_qc(seqs)
        self.assertEqual(qc["internal_stops"], 2)
        self.assertEqual(qc["n_codons"], 2)

    def test_terminal_stop_and_gap(self):
        seqs = {"a": "ATG---", "b": "ATG---"}
        qc = codon_qc(seqs)
        self.assertEqual(qc["terminal_stops"], 0)
        self.assertGreater(qc["gap_fraction"], 0)

    def test_length_must_be_multiple_of_three(self):
        with self.assertRaises(ValueError):
            codon_qc({"a": "ATGA", "b": "ATGA"})


@unittest.skipUnless(importlib.util.find_spec("dendropy"), "dendropy is not installed")
class TestParseAbsrel(unittest.TestCase):
    def test_selected_and_background(self):
        from screen_gene import parse_absrel

        payload = {
            "branch attributes": {
                "0": {
                    "n00001": {
                        "original name": "n00001",
                        "Corrected P-value": 0.01,
                        "Uncorrected P-value": 0.02,
                        "LRT": 5.5,
                        "Baseline MG94xREV omega ratio": 0.4,
                        "Rate Distributions": [[0.2, 0.9], [2.5, 0.1]],
                    },
                    "bat1": {
                        "original name": "bat1",
                        "Corrected P-value": None,
                        "Uncorrected P-value": None,
                        "LRT": 0.0,
                        "Baseline MG94xREV omega ratio": 0.2,
                        "Rate Distributions": [[0.2, 1.0]],
                    },
                }
            },
            "tested": {"0": {"n00001": "test", "bat1": "background"}},
            "analysis": {"version": "2.5.8"},
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "gene.json"
            path.write_text(json.dumps(payload))
            df, version, failures = parse_absrel(
                path, {"n00001": "n00001", "bat1": "bat1"}, "gene", 0.05
            )
        self.assertEqual(version, "2.5.8")
        self.assertEqual(failures, [])
        row = df.set_index("node_label")
        self.assertTrue(bool(row.loc["n00001", "selected"]))
        self.assertAlmostEqual(float(row.loc["n00001", "omega_max"]), 2.5)
        self.assertFalse(bool(row.loc["bat1", "selected"]))
        self.assertFalse(bool(row.loc["bat1", "tested"]))


if __name__ == "__main__":
    unittest.main()
