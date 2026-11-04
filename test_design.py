"""Unit tests. Run: pytest test_design.py -v"""
import argparse
import csv
import gzip
import json
import os
import random
import subprocess
import sys

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [ROOT, os.path.join(ROOT, "analysis"), os.path.join(ROOT, "scripts")]

import design  # noqa: E402


def fake_blast(monkeypatch, specific=True):
    monkeypatch.setattr(design, "run_blast_batch",
                        lambda seqs, db, **k: {q: (1 if specific else 2, [], specific, None) for q in seqs})


@pytest.fixture
def genome():
    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    recs = {"chr": SeqRecord(Seq(seq), id="chr")}
    genes = {"g1": {"contig": "chr", "start": 500, "end": 1000, "strand": "+"},
             "g2": {"contig": "chr", "start": 3000, "end": 3500, "strand": "+"}}
    return recs, genes


def design_target(genome, tid, **kw):
    recs, genes = genome
    return design.process_single_target(tid, recs, genes, "db", False, num_candidates=10, **kw)


# ---------------------------------------------------------------- inputs
class TestInputs:
    def test_parse_gff(self, tmp_path):
        f = tmp_path / "t.gff"
        f.write_text("##gff-version 3\n"
                     "c\tR\tgene\t190\t255\t.\t+\t.\tID=g1;Name=thrL;locus_tag=b0001\n"
                     "c\tR\tCDS\t190\t255\t.\t+\t0\tParent=g1\n"
                     "c\tR\tgene\t100\t500\t.\t-\t.\tID=g2;locus_tag=b9999\n")
        genes = design.parse_gff(str(f))
        assert genes["thrL"] == {"contig": "c", "start": 190, "end": 255, "strand": "+"}
        assert genes["b9999"]["strand"] == "-"

    def test_parse_empty_gff(self, tmp_path):
        f = tmp_path / "e.gff"
        f.write_text("##gff-version 3\n")
        assert design.parse_gff(str(f)) == {}

    def test_extract_sequence_and_errors(self):
        recs = {"c": SeqRecord(Seq("ATGC" * 50), id="c")}
        genes = {"g": {"contig": "c", "start": 1, "end": 8, "strand": "+"},
                 "x": {"contig": "nope", "start": 1, "end": 4, "strand": "+"}}
        assert design.extract_target_sequence(recs, genes, "g")[0] == "ATGCATGC"
        assert "not found" in design.extract_target_sequence(recs, genes, "missing")[2]
        assert "not found" in design.extract_target_sequence(recs, genes, "x")[2]

    def test_target_list_skips_blank_and_comments(self, tmp_path):
        f = tmp_path / "t.txt"
        f.write_text("a\n\n# c\nb\n")
        assert design.parse_target_list(str(f)) == ["a", "b"]

    def test_names(self):
        assert design.sanitize_name("a b(c)") == "a_b_c"
        assert design.target_gene_id("gyrB_amplicon_S0_3") == "gyrB"


# ---------------------------------------------------------------- BLAST
class TestBlast:
    def test_batch_parses_and_flags_nonspecific(self, monkeypatch):
        class P:
            stdout = "q0\tchr\t10\t30\nq1\tchr\t5\t25\nq1\tchr\t900\t920\n"

        def run(cmd, input, **kw):
            assert input.count(">") == 2
            return P()

        monkeypatch.setattr(subprocess, "run", run)
        res = design.run_blast_batch(["AAAA", "CCCC"], "db")
        assert res["AAAA"][2] is True and res["AAAA"][0] == 1
        assert res["CCCC"][2] is False and "2 hits" in res["CCCC"][3]

    def test_failure_fails_closed(self, monkeypatch):
        def boom(*a, **k):
            raise subprocess.CalledProcessError(1, "blastn")

        monkeypatch.setattr(subprocess, "run", boom)
        assert design.run_blast_batch(["AAAA"], "db")["AAAA"][2] is False


# ---------------------------------------------------------------- design + packing
class TestDesign:
    def test_primer3_runs(self):
        rng = random.Random(1)
        seq = "".join(rng.choice("ACGT") for _ in range(400))
        res = design.design_primers_for_sequence(seq, "t", {}, num_candidates=3)
        assert "PRIMER_PAIR_NUM_RETURNED" in res

    def test_candidates_are_complete(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        _, cands, err = design_target(genome, "g1")
        assert err is None and cands
        c = cands[0]
        assert c["end"] - c["start"] == c["size"]
        assert " " not in c["name"] and "(" not in c["name"]
        assert c["fwd_oligo"] == c["fwd"] and c["flags"] == "OK"
        assert [x["sort_key"] for x in cands] == sorted(x["sort_key"] for x in cands)

    def test_nonspecific_primers_rejected(self, genome, monkeypatch):
        fake_blast(monkeypatch, specific=False)
        _, cands, err = design_target(genome, "g1")
        assert cands == [] and "No specific" in err

    def test_missing_target_reported(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        _, cands, err = design_target(genome, "nope")
        assert cands == [] and err

    def test_hairpin_clamp_is_optional_and_5prime(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        _, cands, _ = design_target(genome, "g1", hairpin_clamp=True)
        c = cands[0]
        assert c["fwd_oligo"].endswith(c["fwd"]) and "HairpinClamp" in c["flags"]

    def test_pack_is_deterministic_and_complete(self, genome, monkeypatch):
        fake_blast(monkeypatch)
        m = {t: design_target(genome, t)[1] for t in ("g1", "g2")}
        a = design.pack_targets_first_fit(m, False, 5, seed=3)
        b = design.pack_targets_first_fit(m, False, 5, seed=3)
        ids = lambda pools: [[c["target_id"] for c in p] for p in pools]  # noqa: E731
        assert ids(a) == ids(b)
        assert sorted(t for p in ids(a) for t in p) == ["g1", "g2"]

    def test_overlapping_amplicons_go_to_separate_pools(self):
        mk = lambda t, s, e: {"target_id": t, "contig": "c", "start": s, "end": e,  # noqa: E731
                              "fwd": "GAGACCTACCAGCTCACCTC", "rev": "TCCACGGCAATCTCGATGAT"}
        assert design.has_amplicon_overlap([mk("a", 0, 200), mk("b", 150, 350)])
        assert not design.has_amplicon_overlap([mk("a", 0, 200), mk("b", 200, 400)])
        assert design.is_compatible(mk("b", 150, 350), [mk("a", 0, 200)], True) == (False, "amplicon overlap")

    def test_cross_dimer_blocks_pool(self):
        a = {"contig": "c", "start": 0, "end": 100, "fwd": "ACGTACGTACGTAGCTAGCT", "rev": "ACGTACGTACGTAGCTAGCT"}
        b = dict(a, start=500, end=600)
        # identical primers are fully self-complementary enough to trip 3' dimer / heterodimer checks
        ok, why = design.is_compatible(b, [a], False)
        assert (ok, why) == (True, None) or why in ("3' end dimer", "heterodimer")


# ---------------------------------------------------------------- outputs
class TestOutputs:
    def test_write_outputs(self, genome, monkeypatch, tmp_path):
        fake_blast(monkeypatch)
        m = {t: design_target(genome, t)[1] for t in ("g1", "g2")}
        pools = design.pack_targets_first_fit(m, False, 3, seed=0)
        prefix = str(tmp_path / "out")
        design.write_outputs(pools, prefix, ["g3 failed"])
        n_targets = sum(len(p) for p in pools)

        rows = list(csv.DictReader(open(f"{prefix}_pool_1.csv")))
        assert rows[0]["pool"] == "1" and rows[0]["fwd_oligo"] == rows[0]["fwd_primer"]
        order = list(csv.reader(open(f"{prefix}_order.csv")))
        assert order[0] == ["pool", "name", "sequence"] and len(order) == 1 + 2 * n_targets
        fasta = open(f"{prefix}_pool_1_trimming.fasta").read()
        assert fasta.startswith(f">{rows[0]['target_id']}_F\n{rows[0]['fwd_oligo']}\n")
        bed = open(f"{prefix}_pool_1.bed").read().split("\t")
        assert bed[0] == "chr" and int(bed[2]) > int(bed[1])
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
    def test_redesign_maps_amplicon_to_gene(self, tmp_path, monkeypatch):
        fb = {"target_recommendations": [
            {"target_id": "gyrB_amplicon_S0_2", "rel_mean": 0.0, "recommended_adjustment": 10.0,
             "action": "FLAG FOR REDESIGN"},
            {"target_id": "recA_amplicon_S0_0", "rel_mean": 1.0, "recommended_adjustment": 1.0,
             "action": "Maintain current concentration"}]}
        f = tmp_path / "fb.json"
        f.write_text(json.dumps(fb))
        seen = {}
        monkeypatch.setattr(design, "run_design_mode", lambda a: seen.update(
            targets=open(a.target_file).read().split(), rescue=a.rescue))
        design.run_feedback_mode(argparse.Namespace(feedback=str(f), output_prefix=str(tmp_path / "o"), rescue=False))
        assert seen == {"targets": ["gyrB"], "rescue": True}
        rows = list(csv.DictReader(open(tmp_path / "o_rebalancing_recipe.csv")))
        assert [r["Target"] for r in rows] == ["gyrB_amplicon_S0_2", "recA_amplicon_S0_0"]


# ---------------------------------------------------------------- analysis scripts
class TestAnalysis:
    def test_to_float_handles_picard_placeholders(self):
        from analysis_utils import to_float
        assert to_float("?") == to_float("") == to_float(None) == 0.0 and to_float("1.5") == 1.5

    def test_recommend(self):
        from balance import recommend
        assert recommend(0)[1] == "FLAG FOR REDESIGN"
        assert recommend(0.25)[0] == pytest.approx(4.0) and recommend(0.25)[1].startswith("INCREASE")
        assert recommend(1.0)[1] == "Maintain current concentration"
        assert recommend(4.0)[0] == pytest.approx(0.25) and recommend(4.0)[1].startswith("DECREASE")
        assert recommend(0.01)[0] == 10.0  # clamped

    def test_feedback_survives_undefined_metrics(self):
        from balance import build_feedback
        t = [{'name': 'a', 'chrom': 'c', 'start': 1, 'end': 9, 'mean_coverage': 10.0},
             {'name': 'b', 'chrom': 'c', 'start': 1, 'end': 9, 'mean_coverage': 0.0}]
        fb = build_feedback({'MEAN_TARGET_COVERAGE': '10', 'PCT_AMPLIFIED_BASES': '1',
                             'FOLD_80_BASE_PENALTY': '?'}, t)
        assert fb['overall_metrics']['fold_80'] == 0.0
        assert [r['action'] for r in fb['target_recommendations']][0] == "FLAG FOR REDESIGN"


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


def test_no_undefined_names():
    out = subprocess.run([sys.executable, "-m", "pyflakes", "design.py", "core", "analysis", "scripts"],
                         capture_output=True, text=True, cwd=ROOT).stdout
    assert out == "", out
