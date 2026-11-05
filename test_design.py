"""Unit tests. Run: pytest test_design.py -v"""
import argparse
import csv
import gzip
import json
import os
import random
import shutil
import subprocess
import sys

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [ROOT, os.path.join(ROOT, "analysis"), os.path.join(ROOT, "scripts")]

import design  # noqa: E402


DBS = [("db", "db")]


def fake_blast(monkeypatch, count=1, locs=()):
    monkeypatch.setattr(design, "run_blast_batch",
                        lambda seqs, db, **k: {q: (count, list(locs)) for q in seqs})


@pytest.fixture
def genome():
    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    return {"chr": SeqRecord(Seq(seq), id="chr")}


def make_target(genome, name="g1", start=501, end=1000, **extra):
    row = {"name": name, "reference": "chr", "start": str(start), "end": str(end), **extra}
    t, err = design.resolve_target(row, genome)
    assert err is None, err
    return t


def design_target(genome, name="g1", start=501, end=1000, extra=None, **kw):
    t = make_target(genome, name, start, end, **(extra or {}))
    return design.process_single_target(t, genome, DBS, num_candidates=10, **kw)


# ---------------------------------------------------------------- targets table
class TestTargets:
    def test_read_targets_validates(self, tmp_path):
        f = tmp_path / "t.tsv"
        f.write_text("# c\nname\treference\tstart\tend\na\tchr\t1\t100\n")
        assert design.read_targets(str(f))[0]["name"] == "a"
        f.write_text("name\tbogus\na\t1\n")
        with pytest.raises(ValueError, match="unknown column"):
            design.read_targets(str(f))
        f.write_text("name\nA\nA\n")
        with pytest.raises(ValueError, match="duplicate"):
            design.read_targets(str(f))

    def test_resolve_coordinates_and_sequence(self, genome):
        t = make_target(genome, start=101, end=160)
        assert t["seq"] == str(genome["chr"].seq[100:160]).upper() and t["offset"] == 100
        raw, err = design.resolve_target({"name": "r", "sequence": "acgt" * 40}, genome)
        assert err is None and raw["seq"] == "ACGT" * 40 and raw["contig"] == "r" and not raw["coordinate_based"]

    def test_resolve_errors(self, genome):
        for row, msg in [({"name": "x", "reference": "nope", "start": "1", "end": "9"}, "not found"),
                         ({"name": "x", "reference": "chr", "start": "5990", "end": "7000"}, "outside"),
                         ({"name": "x"}, "either"),
                         ({"name": "x", "sequence": "ACGT" * 30, "junctions": "500"}, "junctions"),
                         ({"name": "x", "sequence": "ACGT" * 30, "span": "0-5"}, "within"),
                         ({"name": "x", "sequence": "ACGT" * 30, "amplicon_min": "300", "amplicon_max": "100"}, "amplicon")]:
            t, err = design.resolve_target(row, genome)
            assert t is None and msg in err, (row, err)

    def test_constraint_conversion(self, genome):
        t = make_target(genome, junctions="141", span="10-19,50", avoid="30-39")
        assert t["junctions"] == [139]            # primer3 wants index of last base before the join
        assert t["span"] == [(9, 10), (49, 1)] and t["avoid"] == [(29, 10)]

    def test_names(self):
        assert design.sanitize_name("a b(c)") == "a_b_c"
        assert design.target_name_from_amplicon("gyrB_amplicon_S0_3") == "gyrB"


# ---------------------------------------------------------------- specificity
class TestSpecificity:
    def test_batch_parses_hits(self, monkeypatch):
        class P:
            stdout = "q0\tchr\t10\t30\nq1\tchr\t5\t25\nq1\tchr\t900\t920\n"

        def run(cmd, input, **kw):
            assert input.count(">") == 2
            return P()

        monkeypatch.setattr(subprocess, "run", run)
        res = design.run_blast_batch(["AAAA", "CCCC"], "db")
        assert res["AAAA"] == (1, [("chr", 10, 30)]) and res["CCCC"][0] == 2

    def test_failure_is_unknown_not_zero(self, monkeypatch):
        def boom(*a, **k):
            raise subprocess.CalledProcessError(1, "blastn")

        monkeypatch.setattr(subprocess, "run", boom)
        assert design.run_blast_batch(["AAAA"], "db")["AAAA"] == (None, [])

    def test_predict_products(self):
        f = [("c", 100, 119)]                       # forward-strand hit
        r = [("c", 349, 330)]                       # reverse-strand hit (sstart > send)
        assert design.predict_products(f, r) == [("c", 250)]
        assert design.predict_products(f, r, max_product=200) == []
        assert design.predict_products(f, [("d", 349, 330)]) == []          # other contig
        assert design.predict_products(f, [("c", 50, 69)]) == []             # same strand
        assert design.predict_products([("c", 500, 519)], [("c", 349, 330)]) == []  # not facing


# ---------------------------------------------------------------- design
class TestDesign:
    def test_primer3_runs(self):
        rng = random.Random(1)
        seq = "".join(rng.choice("ACGT") for _ in range(400))
        assert "PRIMER_PAIR_NUM_RETURNED" in design.design_primers_for_sequence(seq, "t", 3)

    def test_candidates_are_complete(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        _, cands, err = design_target(genome)
        assert err is None and cands
        c = cands[0]
        assert c["end"] - c["start"] == c["size"] and 501 - 1 <= c["start"] < 1000
        assert " " not in c["name"] and c["fwd_oligo"] == c["fwd"] and c["flags"] == "OK"
        assert [x["sort_key"] for x in cands] == sorted(x["sort_key"] for x in cands)
        assert c["hits"] == "F:1 R:1"

    def test_amplicon_range_is_per_target(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        _, cands, _ = design_target(genome, start=501, end=1000, extra={"amplicon_min": "70", "amplicon_max": "110"})
        assert cands and all(c["strategy"] > 0 or 70 <= c["size"] <= 110 for c in cands)
        assert cands[0]["strategy"] == 0 and 70 <= cands[0]["size"] <= 110

    def test_nonspecific_rejected_unless_max_hits_or_multiprime(self, genome, monkeypatch):
        fake_blast(monkeypatch, count=2)
        assert design_target(genome)[1] == []
        assert design_target(genome, extra={"max_hits": "2"})[1]
        assert design_target(genome, force_multiprime=True)[1]

    def test_every_database_is_checked(self, genome, monkeypatch):
        monkeypatch.setattr(design, "run_blast_batch",
                            lambda seqs, db, **k: {q: ((1 if db == "tx" else 3), []) for q in seqs})
        t = make_target(genome)
        two = [("tx", "tx"), ("genome", "genome")]
        assert design.process_single_target(t, genome, two, num_candidates=10)[1] == []
        assert design.process_single_target(t, genome, two[:1], num_candidates=10)[1]

    def test_blast_failure_never_accepted(self, genome, monkeypatch):
        monkeypatch.setattr(design, "run_blast_batch", lambda seqs, db, **k: {q: (None, []) for q in seqs})
        assert design_target(genome, force_multiprime=True)[1] == []

    def test_missing_target_reported(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        t = make_target(genome)
        t["seq"] = ""
        assert design.process_single_target(t, genome, DBS)[1] == []

    def test_junction_constraint_forces_a_spanning_primer(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        # junction before 1-based base 141 of the target; amplicon range suits a 300 bp target
        _, cands, err = design_target(genome, start=1001, end=1300,
                                      extra={"junctions": "141", "amplicon_min": "70", "amplicon_max": "150"})
        assert cands, err
        j = 140  # 0-based index of the first base after the join, within the target
        for c in cands:
            f_lo, f_hi = c["start"] - 1000, c["start"] - 1000 + len(c["fwd"]) - 1
            r_hi = c["end"] - 1000 - 1
            r_lo = r_hi - len(c["rev"]) + 1
            spans = lambda lo, hi: lo <= j - 1 and hi >= j  # noqa: E731
            assert spans(f_lo, f_hi) or spans(r_lo, r_hi), c["name"]

    def test_span_and_avoid_constraints(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        _, cands, err = design_target(genome, start=1001, end=1300,
                                      extra={"span": "140-160", "avoid": "1-40", "amplicon_min": "70", "amplicon_max": "200"})
        assert cands, err
        for c in cands:
            lo, hi = c["start"] - 1001, c["end"] - 1001   # 0-based half-open within the target
            assert lo <= 139 and hi >= 160                  # amplicon contains 140-160
            assert lo >= 40                                 # no primer in the avoided 1-40

    def test_unsatisfiable_constraint_is_reported(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        _, cands, err = design_target(genome, start=1001, end=1100,
                                      extra={"span": "1-100", "amplicon_min": "40", "amplicon_max": "60"})
        assert cands == [] and "constraints" in err

    def test_rescue_only_for_unconstrained_coordinate_targets(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        real = design.design_primers_for_sequence
        t = make_target(genome, start=501, end=1000)
        # Make the original window undesignable; shifted windows are fine.
        monkeypatch.setattr(design, "design_primers_for_sequence",
                            lambda seq, *a, **k: {} if seq == t["seq"] else real(seq, *a, **k))
        assert design.process_single_target(t, genome, DBS, num_candidates=5)[1] == []
        cands = design.process_single_target(t, genome, DBS, num_candidates=5, rescue=True)[1]
        assert cands and cands[0]["flags"].startswith("Rescued")
        assert 500 - 45 <= cands[0]["start"] <= 1000 + 45
        constrained = make_target(genome, start=501, end=1000, avoid="1-5")
        monkeypatch.setattr(design, "design_primers_for_sequence",
                            lambda seq, *a, **k: {} if seq == constrained["seq"] else real(seq, *a, **k))
        assert design.process_single_target(constrained, genome, DBS, num_candidates=5, rescue=True)[1] == []

    def test_pack_is_deterministic_and_complete(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        m = {n: design_target(genome, n, s, s + 499)[1] for n, s in (("g1", 501), ("g2", 3001))}
        a = design.pack_targets_first_fit(m, False, 5, seed=3)
        b = design.pack_targets_first_fit(m, False, 5, seed=3)
        ids = lambda pools: [[c["target_id"] for c in p] for p in pools]  # noqa: E731
        assert ids(a) == ids(b)
        assert sorted(t for p in ids(a) for t in p) == ["g1", "g2"]

    def test_packing_does_not_depend_on_dict_insertion_order(self, monkeypatch):
        # Workers finish in arbitrary order, so the packer must not shuffle in insertion order.
        seen = []
        real_sample = random.Random.sample
        monkeypatch.setattr(random.Random, "sample",
                            lambda self, pop, k: (seen.append(list(pop)), real_sample(self, pop, k))[1])
        monkeypatch.setattr(design, "is_compatible", lambda *a: (False, "x"))  # force many pools: no early exit
        cand = lambda t: [{"target_id": t, "contig": "c", "start": 0, "end": 1,  # noqa: E731
                           "fwd": "ACGTACGTACGTAGCTAGCT", "rev": "TTGGCCAATTGGCCAATTGG"}]
        names = ["t3", "t1", "t2", "t0"]
        for order in (names, list(reversed(names))):
            design.pack_targets_first_fit({n: cand(n) for n in order}, False, 3, seed=1)
        assert seen and all(pop == sorted(pop) for pop in seen)

    def test_cross_target_products_block_a_pool(self, monkeypatch):
        monkeypatch.setitem(design.CFG, "hetero_dimer", -100.0)   # isolate the product rule from dimer checks
        monkeypatch.setitem(design.CFG, "end_3p", -100.0)
        mk = lambda t, s, e, strand="+": {"target_id": t, "contig": "c", "start": s, "end": e, "strand": strand,  # noqa: E731
                                         "fwd": "GAGACCTACCAGCTCACCTC", "rev": "TCCACGGCAATCTCGATGAT"}
        a = mk("a", 0, 200)
        # b is downstream, no overlap: F_a (0) faces R_b (end 700) -> 700 bp product
        assert design.cross_product_size(a, mk("b", 500, 700)) == 700
        monkeypatch.setitem(design.CFG, "max_cross", 1000)
        assert design.is_compatible(mk("b", 500, 700), [a], False)[0] is False
        assert "700 bp" in design.is_compatible(mk("b", 500, 700), [a], False)[1]
        monkeypatch.setitem(design.CFG, "max_cross", 600)
        assert design.is_compatible(mk("b", 500, 700), [a], False)[0] is True
        monkeypatch.setitem(design.CFG, "max_cross", 0)                # disabled
        assert design.is_compatible(mk("b", 500, 700), [a], False)[0] is True
        # a product whose size equals an amplicon length must still be caught (no accidental exclusion)
        monkeypatch.setitem(design.CFG, "max_cross", 1000)
        assert design.is_compatible(mk("b", 100, 200), [mk("a", 0, 200)], False)[0] is False
        # minus-strand targets: roles swap on the reference but the rule is symmetric
        assert design.cross_product_size(mk("a", 0, 200, "-"), mk("b", 500, 700, "-")) is not None
        # another contig never interacts
        other = dict(mk("b", 500, 700), contig="d")
        assert design.is_compatible(other, [a], False)[0] is True

    def test_overlapping_amplicons_go_to_separate_pools(self):
        mk = lambda t, s, e: {"target_id": t, "contig": "c", "start": s, "end": e,  # noqa: E731
                              "fwd": "GAGACCTACCAGCTCACCTC", "rev": "TCCACGGCAATCTCGATGAT"}
        assert design.has_amplicon_overlap([mk("a", 0, 200), mk("b", 150, 350)])
        assert not design.has_amplicon_overlap([mk("a", 0, 200), mk("b", 200, 400)])
        assert design.is_compatible(mk("b", 150, 350), [mk("a", 0, 200)], True) == (False, "amplicon overlap")


# ---------------------------------------------------------------- tiling
class TestTiling:
    F, R = "GAGACCTACCAGCTCACCTC", "TCCACGGCAATCTCGATGAT"

    def cands(self, windows, contig="c", amp=200):
        """windows: {name: (window_start, window_end)}; amplicon = first `amp` bp of each window."""
        return {n: [{"target_id": n, "contig": contig, "start": w[0], "end": w[0] + amp, "strand": "+",
                     "fwd": self.F, "rev": self.R, "tstart": w[0], "tend": w[1]}]
                for n, w in windows.items()}

    def chain(self):   # A,B,C,D in tandem: each window overlaps only its neighbours
        return {"A": (0, 300), "B": (200, 500), "C": (400, 700), "D": (600, 900)}

    @pytest.fixture(autouse=True)
    def only_structural_rules(self, monkeypatch):
        monkeypatch.setitem(design.CFG, "hetero_dimer", -100.0)   # isolate layout rules from dimer checks
        monkeypatch.setitem(design.CFG, "end_3p", -100.0)
        monkeypatch.setitem(design.CFG, "self_dimer", -100.0)

    def test_tandem_targets_alternate_between_two_pools(self):
        assignment, unassigned = design.assign_tiled_pools(self.cands(self.chain()))
        assert assignment == {"A": 0, "B": 1, "C": 0, "D": 1} and unassigned == []

    def test_three_mutually_overlapping_windows_cannot_be_two_coloured(self):
        c = self.cands({"A": (0, 300), "B": (100, 400), "C": (200, 500)})
        assignment, unassigned = design.assign_tiled_pools(c)
        assert unassigned == ["C"] and assignment == {"A": 0, "B": 1}

    def test_isolated_targets_go_to_pool_zero(self):
        c = self.cands({"A": (0, 300), "B": (5000, 5300), "C": (9000, 9300)})
        assert set(design.assign_tiled_pools(c)[0].values()) == {0}

    def test_layout_detection(self):
        assert design.detect_layout(self.cands(self.chain())) == "tiled"
        assert design.detect_layout(self.cands({"A": (0, 300), "B": (5000, 5300)})) == "sparse"
        # windows do not overlap but the chosen amplicons do (amplicon longer than its window slot)
        assert design.detect_layout(self.cands({"A": (0, 150), "B": (160, 310)}, amp=250)) == "tiled"

    def test_fixed_pools_hold_A_C_and_B_D(self, monkeypatch):
        monkeypatch.setitem(design.CFG, "max_cross", 250)           # tiled default = max amplicon size
        c = self.cands(self.chain())
        assignment, _ = design.assign_tiled_pools(c)
        pools, unplaced = design.pack_fixed_pools(c, assignment, iterations=3)
        assert unplaced == [] and [sorted(x["target_id"] for x in p) for p in pools] == [["A", "C"], ["B", "D"]]

    def test_default_cross_product_limit_would_block_tiling(self, monkeypatch):
        # A and C share a pool: A's forward primer faces C's reverse primer ~600 bp away.
        c = self.cands(self.chain())
        assert design.cross_product_size(c["A"][0], c["C"][0]) == 600
        monkeypatch.setitem(design.CFG, "max_cross", 1000)           # the sparse default
        assignment, _ = design.assign_tiled_pools(c)
        _, unplaced = design.pack_fixed_pools(c, assignment, iterations=2)
        assert unplaced                                              # tiling is impossible with it...
        monkeypatch.setitem(design.CFG, "max_cross", 250)            # ...and fine with the tiled default
        assert design.pack_fixed_pools(c, assignment, iterations=2)[1] == []

    def test_tiling_gaps_after_primer_trimming(self):
        def pool(*spans, strand="+"):
            return [[{"contig": "c", "start": s, "end": e, "strand": strand, "fwd": self.F, "rev": self.R}
                     for s, e in spans]]
        # neighbours overlap by 20 bp = one primer length: primers trimmed from both leave a 20 bp hole
        assert design.tiling_gaps(pool((0, 250), (230, 480))) == [("c", 230, 250)]
        # overlap of 40 bp (two primer lengths) leaves no hole; neither does a bigger one
        assert design.tiling_gaps(pool((0, 250), (210, 460))) == []
        assert design.tiling_gaps(pool((0, 250), (150, 400))) == []
        # a real gap between amplicons is reported, and minus-strand amplicons behave the same
        assert design.tiling_gaps(pool((0, 250), (300, 550), strand="-")) == [("c", 230, 320)]
        # the outer ends of the region are not gaps
        assert design.tiling_gaps(pool((0, 250))) == []

    def test_amplicons_overlapping_inside_one_pool_are_not_allowed(self, monkeypatch):
        monkeypatch.setitem(design.CFG, "max_cross", 0)
        c = self.cands({"A": (0, 300), "C": (100, 400)}, amp=250)     # A and C overlap but were assigned together
        _, unplaced = design.pack_fixed_pools(c, {"A": 0, "C": 0}, iterations=2)
        assert unplaced == ["C"]


# ---------------------------------------------------------------- strand and tails
class TestStrandAndTails:
    START, END = 1001, 1300                      # 1-based target region on the reference

    def _cands(self, genome, monkeypatch, strand, **extra):
        fake_blast(monkeypatch)
        extra = dict(extra, strand=strand, amplicon_min="70", amplicon_max="150")
        _, cands, err = design_target(genome, "t", self.START, self.END, extra=extra)
        assert cands, err
        return cands

    @pytest.mark.parametrize("strand", ["+", "-"])
    def test_coordinates_match_sequence_on_either_strand(self, genome, monkeypatch, strand):
        ref = str(genome["chr"].seq).upper()
        for c in self._cands(genome, monkeypatch, strand):
            window = ref[c["start"]:c["end"]]
            sense = window if strand == "+" else design.revcomp(window)
            assert sense.startswith(c["fwd"]), (strand, c["name"])                       # forward = sense primer
            assert sense.endswith(design.revcomp(c["rev"])), (strand, c["name"])          # reverse = antisense primer
            assert c["strand"] == strand and 1000 <= c["start"] and c["end"] <= 1300

    def test_minus_strand_primer_roles_are_swapped_on_the_reference(self, genome, monkeypatch):
        plus = self._cands(genome, monkeypatch, "+")[0]
        minus = self._cands(genome, monkeypatch, "-")[0]
        ref = str(genome["chr"].seq).upper()
        assert ref[plus["start"]:plus["start"] + len(plus["fwd"])] == plus["fwd"]     # fwd leftmost on '+'
        end = minus["end"]
        assert design.revcomp(ref[end - len(minus["fwd"]):end]) == minus["fwd"]       # fwd rightmost on '-'

    def test_minus_strand_constraints_are_in_reference_orientation(self, genome, monkeypatch):
        # (a primer cannot overlap a 'span' range, so the junction sits outside it)
        cands = self._cands(genome, monkeypatch, "-", junctions="201", span="140-160", avoid="1-40")
        a, b = 1000 + 199, 1000 + 200                      # reference indices either side of the join
        for c in cands:
            assert c["start"] <= 1000 + 139 and c["end"] >= 1000 + 160                # amplicon contains 140-160
            f = (c["end"] - len(c["fwd"]), c["end"])                                    # primer footprints on the ref
            r = (c["start"], c["start"] + len(c["rev"]))
            assert all(hi <= 1000 or lo >= 1040 for lo, hi in (f, r))                  # nothing in avoided 1-40
            assert any(lo <= a and hi > b for lo, hi in (f, r)), c["name"]            # a primer spans the join

    def test_bad_strand_rejected(self, genome):
        t, err = design.resolve_target({"name": "x", "sequence": "ACGT" * 40, "strand": "?"}, genome)
        assert t is None and "strand" in err

    def test_tails_are_prepended_and_not_in_trimming_file(self, genome, monkeypatch, tmp_path):
        fake_blast(monkeypatch)
        t = make_target(genome, "t", 1001, 1300)
        _, cands, _ = design.process_single_target(t, genome, DBS, num_candidates=5,
                                                   tails=(design.P5_TAIL, design.P7_TAIL))
        c = cands[0]
        assert c["fwd_oligo"] == design.P5_TAIL + c["fwd"] and c["rev_oligo"] == design.P7_TAIL + c["rev"]
        design.write_outputs([[c]], str(tmp_path / "o"), [])
        order = list(csv.reader(open(tmp_path / "o_order.csv")))
        assert order[1][2] == c["fwd_oligo"] and order[2][2] == c["rev_oligo"]
        fasta = open(tmp_path / "o_pool_1_trimming.fasta").read()
        assert c["fwd"] in fasta and design.P5_TAIL not in fasta       # reads start at the specific primer

    def test_tail_cli_defaults_and_validation(self):
        parser = design.build_parser()
        assert design.resolve_tails(parser.parse_args([])) == ("ACACGACGCTCTTCCGATCT", "GACGTGTGCTCTTCCGATCT")
        assert design.resolve_tails(parser.parse_args(["--no-tails"])) == ("", "")
        assert design.resolve_tails(parser.parse_args(["--fwd-tail", "acgt", "--rev-tail", "ttgg"])) == ("ACGT", "TTGG")
        with pytest.raises(ValueError):
            design.resolve_tails(parser.parse_args(["--fwd-tail", "ACGN"]))


# ---------------------------------------------------------------- outputs
class TestOutputs:
    def test_write_outputs(self, genome, monkeypatch, tmp_path):
        fake_blast(monkeypatch)
        m = {n: design_target(genome, n, s, s + 499)[1] for n, s in (("g1", 501), ("g2", 3001))}
        pools = design.pack_targets_first_fit(m, False, 3, seed=0)
        prefix = str(tmp_path / "out")
        design.write_outputs(pools, prefix, ["g3 failed"])
        n_targets = sum(len(p) for p in pools)

        rows = list(csv.DictReader(open(f"{prefix}_pool_1.csv")))
        assert rows[0]["pool"] == "1" and rows[0]["fwd_oligo"] == rows[0]["fwd_primer"]
        order = list(csv.reader(open(f"{prefix}_order.csv")))
        assert order[0] == ["pool", "name", "sequence", "role"] and len(order) == 1 + 2 * n_targets
        assert order[1][3].startswith("forward") and "RT primer" in order[2][3]
        fasta = open(f"{prefix}_pool_1_trimming.fasta").read()
        assert fasta.startswith(f">{rows[0]['target_id']}_F\n{rows[0]['fwd_oligo']}\n")
        bed = open(f"{prefix}_pool_1.bed").read().split("\t")
        assert bed[0] == "chr" and int(bed[2]) > int(bed[1])
        assert "predicted_products" in rows[0]
        assert "g3 failed" in open(f"{prefix}_failed_targets.txt").read()


# ---------------------------------------------------------------- biophysics
class TestBiophysics:
    def test_invalid_sequence_not_reported_as_no_dimer(self):
        from core.biophysics import get_homodimer_dg
        assert get_homodimer_dg("ACGT" * 40) < -100  # >60 nt is rejected by primer3

    def test_heterodimer_symmetric(self):
        from core.biophysics import get_heterodimer_dg
        a, b = "GAGACCTACCAGCTCACCTC", "TCCACGGCAATCTCGATGAT"
        assert get_heterodimer_dg(a, b) == get_heterodimer_dg(b, a)


# ---------------------------------------------------------------- feedback loop
class TestFeedback:
    FB = {"target_recommendations": [
        {"target_id": "gyrB_amplicon_S0_2", "rel_mean": 0.0, "recommended_adjustment": 10.0,
         "action": "FLAG FOR REDESIGN"},
        {"target_id": "recA_amplicon_S0_0", "rel_mean": 1.0, "recommended_adjustment": 1.0,
         "action": "Maintain current concentration"}]}

    def _args(self, tmp_path, targets):
        (tmp_path / "fb.json").write_text(json.dumps(self.FB))
        return argparse.Namespace(feedback=str(tmp_path / "fb.json"), output_prefix=str(tmp_path / "o"),
                                  rescue=False, targets=targets)

    def test_redesign_selects_flagged_target_rows(self, tmp_path, monkeypatch):
        tsv = tmp_path / "panel.tsv"
        tsv.write_text("name\treference\tstart\tend\ngyrB\tchr\t1\t500\nrecA\tchr\t600\t900\n")
        seen = {}
        monkeypatch.setattr(design, "run_design_mode", lambda a: seen.update(
            rows=design.read_targets(a.targets), rescue=a.rescue, prefix=a.output_prefix))
        design.run_feedback_mode(self._args(tmp_path, str(tsv)))
        assert [r["name"] for r in seen["rows"]] == ["gyrB"] and seen["rescue"] is True
        assert seen["prefix"].endswith("o_redesign")
        recipe = list(csv.DictReader(open(tmp_path / "o_rebalancing_recipe.csv")))
        assert [r["Target"] for r in recipe] == ["gyrB_amplicon_S0_2", "recA_amplicon_S0_0"]

    def test_without_targets_only_recipe_is_written(self, tmp_path, monkeypatch):
        monkeypatch.setattr(design, "run_design_mode", lambda a: pytest.fail("must not redesign"))
        design.run_feedback_mode(self._args(tmp_path, None))
        assert (tmp_path / "o_rebalancing_recipe.csv").exists()

    def test_names_with_odd_characters_map_back(self, tmp_path, monkeypatch):
        tsv = tmp_path / "panel.tsv"
        tsv.write_text("name\treference\tstart\tend\ngyrB (v2)\tchr\t1\t500\n")
        fb = {"target_recommendations": [{"target_id": "gyrB_v2_amplicon_S0_1", "rel_mean": 0.0,
                                          "recommended_adjustment": 10.0, "action": "FLAG FOR REDESIGN"}]}
        (tmp_path / "fb.json").write_text(json.dumps(fb))
        seen = {}
        monkeypatch.setattr(design, "run_design_mode", lambda a: seen.update(rows=design.read_targets(a.targets)))
        design.run_feedback_mode(argparse.Namespace(feedback=str(tmp_path / "fb.json"), targets=str(tsv),
                                                    output_prefix=str(tmp_path / "o"), rescue=False))
        assert seen["rows"][0]["name"] == "gyrB (v2)"


# ---------------------------------------------------------------- analysis (samtools + Python)
class TestAnalysis:
    def test_recommend(self):
        from balance import recommend
        assert recommend(0)[1] == "FLAG FOR REDESIGN"
        assert recommend(0.25)[0] == pytest.approx(4.0) and recommend(0.25)[1].startswith("INCREASE")
        assert recommend(1.0)[1] == "Maintain current concentration"
        assert recommend(4.0)[0] == pytest.approx(0.25) and recommend(4.0)[1].startswith("DECREASE")
        assert recommend(0.01)[0] == 10.0  # clamped

    def test_feedback_handles_undefined_fold80_and_flags_zero_coverage(self):
        from balance import build_feedback
        t = [{"name": "a", "chrom": "c", "start": 1, "end": 9, "mean_coverage": 10.0},
             {"name": "b", "chrom": "c", "start": 1, "end": 9, "mean_coverage": 0.0}]
        fb = build_feedback({"mean_target_coverage": 10.0, "pct_on_target_bases": 99.0, "fold_80": None}, t)
        assert fb["overall_metrics"]["fold_80"] is None and fb["overall_metrics"]["uniformity_0.2x"] == 50.0
        assert fb["target_recommendations"][0]["action"] == "FLAG FOR REDESIGN"
        json.dumps(fb)                                                       # serialisable (None -> null)

    def test_parse_bedcov_and_stats(self):
        import metrics
        six = "chr\t100\t200\tT1_amplicon_S0_0\t0\t+\t5000\nchr\t300\t350\tT2\t0\t-\t0\n"
        t = metrics.parse_bedcov(six)
        assert [(x["name"], x["mean_coverage"]) for x in t] == [("T1_amplicon_S0_0", 50.0), ("T2", 0.0)]
        assert metrics.parse_bedcov("chr\t0\t10\t40\n")[0]["name"] == "chr:0-10"      # 3-column BED gets a name
        stats = ("SN\traw total sequences:\t200\t# x\nSN\treads mapped:\t190\t# x\n"
                 "SN\tbases mapped (cigar):\t17000\t# x\nSN\treads paired:\t200\t# x\n")
        assert metrics.parse_stats(stats) == {"reads_total": 200, "reads_mapped": 190, "bases_mapped": 17000}

    def test_fold_80_matches_the_definition(self):
        import numpy as np
        import metrics
        d = np.array([100] * 80 + [50] * 20)
        assert metrics.fold_80(d) == pytest.approx(d.mean() / np.percentile(d, 20))
        assert metrics.fold_80([5, 5, 5, 5]) == pytest.approx(1.0)
        assert metrics.fold_80([0, 0, 0, 10]) is None and metrics.fold_80([]) is None   # undefined, Picard prints '?'

    def test_panel_metrics(self):
        import numpy as np
        import metrics
        m = metrics.panel_metrics([{}, {}], np.array([100.0] * 50 + [300.0] * 50),
                                  {"reads_total": 100, "reads_mapped": 90, "bases_mapped": 25000})
        assert m["mean_target_coverage"] == 200.0 and m["target_territory_bp"] == 100
        assert m["pct_on_target_bases"] == pytest.approx(80.0) and m["pct_reads_mapped"] == 90.0

    def test_find_samples(self, tmp_path):
        import analyze
        for n in ("A_R1.fastq.gz", "A_R2.fastq.gz", "B_R1.fastq.gz", "B_R2.fastq.gz"):
            (tmp_path / n).write_text("")
        got = analyze.find_samples(str(tmp_path / "*_R{1,2}.fastq.gz"))
        assert [s for s, _ in got] == ["A", "B"] and all(len(f) == 2 for _, f in got)
        assert got[0][1][1].endswith("A_R2.fastq.gz")
        se = analyze.find_samples(str(tmp_path / "A_R1.fastq.gz"))
        assert se == [("A_R1", [str(tmp_path / "A_R1.fastq.gz")])]
        (tmp_path / "C_R1.fastq.gz").write_text("")
        with pytest.raises(FileNotFoundError, match="mate"):
            analyze.find_samples(str(tmp_path / "*_R{1,2}.fastq.gz"))
        with pytest.raises(ValueError, match="wildcards"):
            analyze.find_samples(str(tmp_path / "A_R{1,2}*.fastq.gz"))

    def test_split_primers(self, tmp_path):
        import analyze
        fa = tmp_path / "p.fa"
        fa.write_text(">a_F\nACGT\n>a_R\nTTTT\n>b_F\nGGGG\n>b_R\nCCCC\n")
        f, r = analyze.split_primers(str(fa), str(tmp_path))
        assert open(f).read() == ">a_F\nACGT\n>b_F\nGGGG\n" and open(r).read() == ">a_R\nTTTT\n>b_R\nCCCC\n"
        fa.write_text(">x\nACGT\n")
        with pytest.raises(ValueError, match="no forward primers"):
            analyze.split_primers(str(fa), str(tmp_path))

    @pytest.mark.skipif(not all(shutil.which(t) for t in ("cutadapt", "bwa", "samtools")),
                        reason="needs cutadapt, bwa and samtools (the conda env)")
    def test_end_to_end_on_the_small_fixture(self, tmp_path):
        import analyze
        import make_test_fixture as fx
        fx.build(os.path.join(ROOT, "common_refs", "ecoli_genome.fna"), str(tmp_path / "fx"), coverage=40)
        fx_dir = str(tmp_path / "fx")
        args = (os.path.join(fx_dir, "FixtureSample_R{1,2}.fastq.gz"), os.path.join(fx_dir, "genome.fna"),
                os.path.join(fx_dir, "targets.bed"), os.path.join(fx_dir, "primers_trimming.fasta"), str(tmp_path / "res"))
        (s,) = analyze.analyze(*args, threads=2)
        assert s["sample"] == "FixtureSample" and s["n_targets"] == 3
        assert s["pct_retained"] == 100.0 and s["pct_reads_mapped"] > 99 and s["pct_on_target_bases"] > 95
        assert s["mean_target_coverage"] > 20
        prov = s["provenance"]                                   # what a re-run needs to reproduce or compare
        assert {"cutadapt", "samtools", "bwa", "python"} <= set(prov["tool_versions"]) and "unknown" not in prov["tool_versions"].values()
        assert all(prov["inputs"][k]["sha256"] for k in ("genome", "bed", "primers")) and len(prov["reads"]) == 2
        for f in ("FixtureSample.sorted.bam", "FixtureSample_balancing_feedback.json", "FixtureSample_target_summary.tsv",
                  "FixtureSample_uniformity.png", "target_coverage.tsv", "cutadapt.log", "FixtureSample_summary.json"):
            assert os.path.exists(os.path.join(tmp_path, "res", "FixtureSample", f)), f
        fb = json.load(open(tmp_path / "res/FixtureSample/FixtureSample_balancing_feedback.json"))
        assert len(fb["target_recommendations"]) == 3
        # a primer FASTA from another pool must fail loudly, not produce an empty report
        bad = tmp_path / "other.fa"
        bad.write_text(">x_F\nACGTACGTACGTACGTACGT\n>x_R\nTTTTGGGGAAAACCCCTTTT\n")
        with pytest.raises(RuntimeError, match="no reads started with a primer"):
            analyze.analyze(args[0], args[1], args[2], str(bad), str(tmp_path / "res2"), threads=2)

    @pytest.mark.skipif(not all(shutil.which(t) for t in ("cutadapt", "bwa", "samtools")),
                        reason="needs cutadapt, bwa and samtools (the conda env)")
    def test_either_orientation_recovers_reads_that_start_at_the_reverse_primer(self, tmp_path):
        import analyze
        import make_test_fixture as fx
        fx.build(os.path.join(ROOT, "common_refs", "ecoli_genome.fna"), str(tmp_path / "fx"), coverage=40)
        d = str(tmp_path / "fx")
        raw = lambda n: open(os.path.join(d, f"FixtureSample_R{n}.fastq.gz"), "rb").read()
        mixed = tmp_path / "mixed"                                # half the pairs as made, half with R1/R2 swapped
        mixed.mkdir()
        (mixed / "M_R1.fastq.gz").write_bytes(raw(1) + raw(2))
        (mixed / "M_R2.fastq.gz").write_bytes(raw(2) + raw(1))
        args = (str(mixed / "M_R{1,2}.fastq.gz"), os.path.join(d, "genome.fna"), os.path.join(d, "targets.bed"),
                os.path.join(d, "primers_trimming.fasta"))
        (default,) = analyze.analyze(*args, str(tmp_path / "a"), threads=2)
        (either,) = analyze.analyze(*args, str(tmp_path / "b"), threads=2, either_orientation=True)
        assert default["pct_retained"] == pytest.approx(50.0, abs=0.5)
        assert either["pct_retained"] == 100.0 and either["pct_reads_mapped"] > 99
        assert either["n_targets"] == 3 and either["mean_target_coverage"] > 1.8 * default["mean_target_coverage"]


class TestSimulatorAndFixture:
    def test_simulator_deterministic_and_writes_truth(self, tmp_path):
        import generate_simulated_fastq as sim
        g = {"c": SeqRecord(Seq("ACGT" * 500), id="c")}
        outs = []
        for n in ("a", "b"):
            sim.generate_fastq_records([("c", 100, 400, "t1")], g, 20, 100, str(tmp_path / n),
                                       seed=5, non_uniform=True)
            outs.append(gzip.open(f"{tmp_path / n}_R1.fastq.gz", "rt").read())
        assert outs[0] == outs[1] and outs[0]
        assert (tmp_path / "a.truth.tsv").read_text().startswith("chrom\tpos\tref\talt\tamplicon")

    def test_fixture_builds_primer_only_fasta(self, tmp_path):
        import make_test_fixture as fx
        out = fx.build(os.path.join(ROOT, "common_refs", "ecoli_genome.fna"), str(tmp_path), coverage=5)
        fasta = open(os.path.join(out, "primers_trimming.fasta")).read().split("\n")
        assert fasta[0] == ">t1_F" and len(fasta[1]) == fx.PRIMER_LEN
        for f in ("genome.fna", "targets.bed", "FixtureSample_R1.fastq.gz", "FixtureSample_R2.fastq.gz"):
            assert os.path.exists(os.path.join(out, f))


class TestProjectAndManifest:
    def write(self, tmp_path, text):
        f = tmp_path / "panel.toml"
        f.write_text(text)
        return str(f)

    def test_project_sets_defaults_and_cli_overrides(self, tmp_path):
        proj = self.write(tmp_path, 'name = "demo"\ntargets = "t.tsv"\nreference = ["a.fa", "b.fa"]\n'
                          'output_prefix = "out/p"\n[design]\nlayout = "tiled"\namplicon_max = 300\nseed = 7\n'
                          'rescue = true\n[tails]\nenabled = false\n[thresholds]\nself_dimer_dg = -10\n')
        _, args = design.parse_args(["--project", proj, "--seed", "9"])
        assert args.targets == str(tmp_path / "t.tsv") and args.output_prefix == str(tmp_path / "out/p")
        assert args.reference == [str(tmp_path / "a.fa"), str(tmp_path / "b.fa")]
        assert (args.layout, args.amplicon_max, args.rescue, args.project_name) == ("tiled", 300, True, "demo")
        assert args.seed == 9                      # CLI wins over the project file
        assert args.no_tails is True and design.resolve_tails(args) == ("", "")
        assert args.max_self_dimer_dg == -10

    def test_no_project_keeps_defaults(self):
        _, args = design.parse_args([])
        assert args.layout == "auto" and args.no_tails is False and args.max_cross_product is None

    @pytest.mark.parametrize("text,msg", [
        ('bogus = 1\n', "unknown key 'bogus'"),
        ('[nope]\nx = 1\n', "unknown section"),
        ('[design]\nseed = "x"\n', "type int"),
        ('[design]\nrescue = 1\n', "type bool"),
        ('reference = "a.fa"\n', "type paths"),
        ('[design]\nmax_cross = 5\n', "unknown key 'max_cross'")])
    def test_bad_project_is_rejected_with_a_clear_message(self, tmp_path, text, msg):
        with pytest.raises(ValueError, match=msg):
            design.load_project(self.write(tmp_path, text))

    def test_example_project_file_is_valid_and_points_at_real_files(self):
        _, args = design.parse_args(["--project", os.path.join(ROOT, "panel.example.toml")])
        assert os.path.exists(args.targets) and all(os.path.exists(r) for r in args.reference)
        assert args.layout == "auto" and design.resolve_tails(args) == (design.P5_TAIL, design.P7_TAIL)

    def test_file_info_hashes_small_files_only(self, tmp_path):
        import analyze
        f = tmp_path / "x"
        f.write_bytes(b"abc")
        assert analyze.file_info(str(f))["sha256"] == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        assert analyze.file_info(str(f), hash_limit=2)["sha256"] is None

    def test_manifest_records_inputs_results_and_next_steps(self, tmp_path, genome, monkeypatch):
        fake_blast(monkeypatch)
        (tmp_path / "t.tsv").write_text("name\treference\tstart\tend\ng1\tchr\t501\t1000\n")
        (tmp_path / "r.fa").write_text(">chr\nACGT\n")
        _, args = design.parse_args(["--targets", str(tmp_path / "t.tsv"), "--reference", str(tmp_path / "r.fa"),
                                     "--output-prefix", str(tmp_path / "o")])
        c = design_target(genome)[1]
        m = design.build_manifest(args, "sparse", [[dict(c[0], rank_used=0)]], ["x failed"], 2,
                                  [("chr", 10, 30)])
        assert m["layout"] == "sparse" and m["results"]["targets_requested"] == 2
        assert m["git_dirty"] in (True, False, None)
        assert {"python", "primer3-py", "biopython", "blast"} <= set(m["tool_versions"])
        assert m["results"]["targets_designed"] == 1 and m["results"]["failed"] == ["x failed"]
        assert m["results"]["pools"] == [{"pool": 1, "n_targets": 1, "targets": ["g1"]}]
        assert m["results"]["tiling_gaps"] == [{"contig": "chr", "start": 10, "end": 30, "length": 20}]
        assert m["inputs"]["targets"]["sha256"] and m["inputs"]["reference"][0]["size_bytes"] == 10
        assert any(o.endswith("_order.csv") for o in m["outputs"]) and any("design_viewer" in n for n in m["next_steps"])
        assert m["options"]["seed"] == 0 and json.dumps(m)       # JSON-serialisable


class TestViewer:
    def _row(self, tid, pool, start, end, strand="+"):
        r = {"target_id": tid, "name": f"{tid}_amplicon_S0_0", "pool": pool, "contig": "c", "start": start,
             "end": end, "strand": strand, "fwd_primer": "A" * 20, "rev_primer": "C" * 20}
        fl, rl = 20, 20
        r["f"], r["r"] = ((start, start + fl), (end - rl, end)) if strand == "+" else ((end - fl, end), (start, start + rl))
        return r

    def test_near_matches_filter_on_3prime_end_and_mismatches(self, monkeypatch):
        import design_viewer as dv
        # qseqid sseqid length mismatch gaps qstart qend sstart send qseq sseq
        lines = "\n".join([
            "p\tc\t17\t1\t0\t4\t20\t100\t116\tAAAAAAAAAAAAACGT\tAAAAAAAAAAAAACGT",   # ok (3' end aligned, last 4 match)
            "p\tc\t17\t1\t0\t1\t17\t300\t316\tAAAAAAAAAAAAACGT\tAAAAAAAAAAAAACGT",   # 3' end not aligned (qend 17 != 20)
            "p\tc\t17\t1\t0\t4\t20\t500\t516\tAAAAAAAAAAAAACGT\tAAAAAAAAAAAAACGA",   # last base mismatch
            "p\tc\t17\t3\t0\t4\t20\t700\t716\tAAAAAAAAAAAAACGT\tAAAAAAAAAAAAACGT",   # too many mismatches
            "p\tc\t12\t0\t0\t9\t20\t900\t911\tAAAAAAAAAAAAACGT\tAAAAAAAAAAAAACGT",   # too short
            "p\tc\t17\t1\t0\t4\t20\t1016\t1000\tAAAAAAAAAAAAACGT\tAAAAAAAAAAAAACGT",  # minus strand, ok
        ]) + "\n"
        monkeypatch.setattr(subprocess, "run", lambda *a, **k: type("R", (), {"stdout": lines})())
        sites = dv.near_matches({"p": "A" * 20}, "db", min_len=16, max_mismatch=2)
        assert [(s["lo"], s["strand"]) for s in sites] == [(99, "+"), (999, "-")]

    def test_own_site_is_excluded_and_unintended_products_are_found(self):
        import design_viewer as dv
        rows = [self._row("a", 1, 0, 200), self._row("b", 1, 5000, 5200), self._row("c", 2, 5000, 5200)]
        sites = [{"primer": "a_F", "contig": "c", "lo": 0, "hi": 20, "strand": "+"},          # its own site
                 {"primer": "b_R", "contig": "c", "lo": 600, "hi": 620, "strand": "-", "len": 17, "mm": 0}]
        mis, products = dv.build(rows, sites, max_product=1000)
        assert [m["primer"] for m in mis] == ["b_R"] and mis[0]["target"] == "b" and mis[0]["pool"] == 1
        # a_F (+, 0-20) faces b_R's near-match site (-, 600-620): a 620 bp product in pool 1
        assert any(p["a"] == "a_F" and p["b"] == "b_R" and p["size"] == 620 and p["pool"] == 1 for p in products)
        # a target's own F+R pair is never reported, and pool 2's lone target has no partner
        assert not any({p["a"], p["b"]} == {"a_F", "a_R"} for p in products)
        assert not any(p["pool"] == 2 for p in products)

    def test_minus_strand_footprints_and_tracks(self, tmp_path):
        import design_viewer as dv
        r = self._row("m", 1, 100, 300, "-")
        assert r["f"] == (280, 300) and r["r"] == (100, 120)
        dv.write_tracks(str(tmp_path), [r], [], [], "ref.fa")
        bed = open(tmp_path / "amplicons.bed").read().splitlines()[1].split("\t")
        assert bed[5] == "-" and bed[9] == "2" and bed[10] == "20,20" and bed[11] == "0,180"
        assert (tmp_path / "session.xml").exists()


class TestEcoliTargets:
    def test_spread_sets_are_genome_wide_unique_and_committed(self):
        import make_ecoli_targets as mk
        genes = mk.read_genes(os.path.join(ROOT, "common_refs", "genomic.gff"), 600)
        assert len(genes) > 2500 and not any(mk.REPEATS.match(g[0]) for g in genes)
        picks = mk.spread(genes, 100)
        starts = [g[2] for g in picks]
        assert len({g[0] for g in picks}) == 100 and starts == sorted(starts)
        assert starts[-1] - starts[0] > 4_000_000 and min(b - a for a, b in zip(starts, starts[1:])) > 20_000
        for n in (50, 100, 250):                        # the committed files are exactly what the generator makes
            rows = [ln.rstrip("\n").split("\t") for ln in open(os.path.join(ROOT, "design", "targets", f"ecoli_spread_{n}.tsv"))][1:]
            assert [r[0] for r in rows] == [g[0] for g in mk.spread(genes, n)]
        with pytest.raises(ValueError, match="eligible"):
            mk.spread(genes, 10**6)

    def test_tiles_overlap_by_window_minus_step(self):
        import make_ecoli_targets as mk
        t = mk.tiled(4, 1000, 300, 200)
        assert [(x[2], x[3]) for x in t] == [(1000, 1299), (1200, 1499), (1400, 1699), (1600, 1899)]
        assert all(x[4:] == (230, 300) for x in t)


def test_genes_to_targets(tmp_path):
    gff = tmp_path / "a.gff"
    gff.write_text("##gff-version 3\nc\tR\tgene\t10\t99\t.\t+\t.\tID=g;Name=geneA;locus_tag=b1\n")
    genes = tmp_path / "g.txt"
    genes.write_text("geneA\n# skip\nmissing\n")
    out = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "genes_to_targets.py"),
                          "--gff", str(gff), "--genes", str(genes)], capture_output=True, text=True)
    assert out.stdout.splitlines() == ["name\treference\tstart\tend", "geneA\tc\t10\t99"]
    assert "missing" in out.stderr


def test_no_undefined_names():
    out = subprocess.run([sys.executable, "-m", "pyflakes", "design.py", "core", "analysis", "scripts"],
                         capture_output=True, text=True, cwd=ROOT).stdout
    assert out == "", out
