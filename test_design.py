"""
Unit tests for design.py core functions.
Run with: pytest test_design.py -v
"""
import pytest
import tempfile
import os
import sys
import importlib.util
from io import StringIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Import the module with dots in filename using importlib
spec = importlib.util.spec_from_file_location("design", "design.py")
design = importlib.util.module_from_spec(spec)
sys.modules["design"] = design
spec.loader.exec_module(design)


class TestParseGff:
    """Tests for the parse_gff function."""
    
    def test_parse_valid_gff(self, tmp_path):
        """Test parsing a valid GFF file with gene entries."""
        gff_content = """##gff-version 3
NC_000913.3	RefSeq	gene	190	255	.	+	.	ID=gene-b0001;Name=thrL;locus_tag=b0001
NC_000913.3	RefSeq	gene	337	2799	.	+	.	ID=gene-b0002;Name=thrA;locus_tag=b0002
NC_000913.3	RefSeq	CDS	337	2799	.	+	0	ID=cds-b0002;Parent=gene-b0002
"""
        gff_file = tmp_path / "test.gff"
        gff_file.write_text(gff_content)
        
        result = design.parse_gff(str(gff_file))
        
        assert "thrL" in result
        assert "thrA" in result
        assert result["thrL"]["start"] == 190
        assert result["thrL"]["end"] == 255
        assert result["thrL"]["strand"] == "+"
        assert result["thrL"]["contig"] == "NC_000913.3"
    
    def test_parse_empty_gff(self, tmp_path):
        """Test parsing an empty GFF file returns empty dict."""
        gff_file = tmp_path / "empty.gff"
        gff_file.write_text("##gff-version 3\n")
        
        result = design.parse_gff(str(gff_file))
        
        assert result == {}
    
    def test_parse_gff_with_locus_tag(self, tmp_path):
        """Test parsing GFF that uses locus_tag instead of Name."""
        gff_content = """##gff-version 3
NC_000913.3	RefSeq	gene	100	500	.	-	.	ID=gene-test;locus_tag=b9999
"""
        gff_file = tmp_path / "locus.gff"
        gff_file.write_text(gff_content)
        
        result = design.parse_gff(str(gff_file))
        
        assert "b9999" in result
        assert result["b9999"]["strand"] == "-"


class TestExtractTargetSequence:
    """Tests for the extract_target_sequence function."""
    
    def test_extract_valid_sequence(self):
        """Test extracting a sequence for a valid target."""
        genome_records = {
            "NC_000913.3": SeqRecord(Seq("ATGCATGCATGCATGCATGC" * 10), id="NC_000913.3")
        }
        gene_coords = {
            "testGene": {"contig": "NC_000913.3", "start": 1, "end": 20, "strand": "+"}
        }
        
        seq, info, error = design.extract_target_sequence(genome_records, gene_coords, "testGene")
        
        assert seq == "ATGCATGCATGCATGCATGC"
        assert info["start"] == 1
        assert error is None
    
    def test_extract_missing_target(self):
        """Test error handling for missing target ID."""
        genome_records = {"NC_000913.3": SeqRecord(Seq("ATGC"), id="NC_000913.3")}
        gene_coords = {}
        
        seq, info, error = design.extract_target_sequence(genome_records, gene_coords, "missingGene")
        
        assert seq is None
        assert info is None
        assert "not found" in error.lower()
    
    def test_extract_missing_contig(self):
        """Test error handling for missing contig."""
        genome_records = {"chr1": SeqRecord(Seq("ATGC"), id="chr1")}
        gene_coords = {
            "testGene": {"contig": "chr2", "start": 1, "end": 4, "strand": "+"}
        }
        
        seq, info, error = design.extract_target_sequence(genome_records, gene_coords, "testGene")
        
        assert seq is None
        assert "not found" in error.lower()


class TestAddIterativeHairpinClamp:
    """Tests for the add_iterative_hairpin_clamp function."""
    
    def test_adds_clamp_to_oligo(self):
        """Test that a clamp sequence is prepended."""
        oligo = "ACGTACGTACGT"
        primer_seq = "GCTAGCTAGCTA"
        target_dg = -10.0
        max_len = 10
        
        result = design.add_iterative_hairpin_clamp(oligo, primer_seq, target_dg, max_len)
        
        # Result should be longer than original (clamp added)
        assert len(result) >= len(oligo)
        # Original sequence should be at the end
        assert result.endswith(oligo)
    
    def test_returns_original_on_short_primer(self):
        """Test returns original oligo if primer is too short."""
        oligo = "ACGT"
        primer_seq = "A"  # Very short
        
        result = design.add_iterative_hairpin_clamp(oligo, primer_seq, -10.0, 10)
        
        # Should still contain the original
        assert oligo in result


class TestBlastDbCreation:
    """Tests for BLAST database creation logic."""
    
    def test_cleaned_fasta_path_no_double_extension(self):
        """Test that .cleaned extension is not doubled."""
        # This tests the logic we fixed earlier
        original_path = "/path/to/genome.fna"
        expected = "/path/to/genome.cleaned.fna"
        
        # Simulate the logic from create_blast_db_if_needed
        if original_path.endswith(".cleaned.fna"):
            cleaned_path = original_path
        else:
            cleaned_path = original_path.replace(".fna", ".cleaned.fna")
        
        assert cleaned_path == expected
        
    def test_already_cleaned_path_unchanged(self):
        """Test that already cleaned path stays the same."""
        original_path = "/path/to/genome.cleaned.fna"
        
        if original_path.endswith(".cleaned.fna"):
            cleaned_path = original_path
        else:
            cleaned_path = original_path.replace(".fna", ".cleaned.fna")
        
        assert cleaned_path == original_path


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


class TestRegressionUndefinedNames:
    """Guards against the missing-import/constant bugs that broke design.py."""

    def test_primer3_design_runs(self):
        import random
        random.seed(1)
        seq = "".join(random.choice("ACGT") for _ in range(400))
        res = design.design_primers_for_sequence(seq, "t", {}, num_candidates=3)
        assert "PRIMER_PAIR_NUM_RETURNED" in res

    def test_adapter_and_hairpin_constants(self):
        for name in ("FWD_P5_TRUNCATED", "REV_P7_TRUNCATED", "FWD_RC_P5_TRUNCATED",
                     "REV_RC_P7_TRUNCATED", "HAIRPIN_STEM_TARGET_DG", "HAIRPIN_CLAMP_MAX_LEN"):
            assert hasattr(design, name)

    def test_write_outputs(self, tmp_path):
        row = {'target_id': 'g', 'pair_rank': '0 (S)', 'flags': 'OK',
               'fwd_primer_tailed': 'A', 'rev_primer_tailed': 'C',
               'fwd_primer_seq': 'ACGTACGTAC', 'rev_primer_seq': 'TTGGCCAATT',
               'fwd_primer_tm': '60', 'rev_primer_tm': '60',
               'amplicon_size': 150, 'specificity_hits': 'F:1, R:1'}
        prefix = str(tmp_path / "out")
        design.write_design_output_files([row], ["c\t0\t150\tg\t0\t+\n"], prefix, [], [])
        assert (tmp_path / "out.csv").exists()
        fasta = (tmp_path / "out_trimming.fasta").read_text()
        assert fasta.startswith(">g_F\n" + design.FWD_P5_TRUNCATED)

    def test_no_undefined_names(self):
        import subprocess
        out = subprocess.run([sys.executable, "-m", "pyflakes", "design.py"],
                             capture_output=True, text=True).stdout
        assert "undefined name" not in out, out


class TestBlastBatch:
    def test_batch_parses_hits_and_flags_nonspecific(self, monkeypatch):
        import subprocess

        class P:
            stdout = "q0\tchr\t10\t30\nq1\tchr\t5\t25\nq1\tchr\t900\t920\n"

        def fake_run(cmd, input, **kw):
            assert input.count(">") == 2
            return P()

        monkeypatch.setattr(subprocess, "run", fake_run)
        res = design.run_blast_batch(["AAAA", "CCCC"], "db")
        assert res["AAAA"][2] is True and res["AAAA"][0] == 1
        assert res["CCCC"][2] is False and "2 hits" in res["CCCC"][3]

    def test_blast_failure_fails_closed(self, monkeypatch):
        import subprocess

        def boom(*a, **k):
            raise subprocess.CalledProcessError(1, "blastn")

        monkeypatch.setattr(subprocess, "run", boom)
        res = design.run_blast_batch(["AAAA"], "db")
        assert res["AAAA"][2] is False


class TestDesignAndPacking:
    @pytest.fixture
    def genome(self):
        import random
        rng = random.Random(7)
        seq = "".join(rng.choice("ACGT") for _ in range(6000))
        return {"chr": SeqRecord(Seq(seq), id="chr")}, {
            "g1": {"contig": "chr", "start": 500, "end": 1000, "strand": "+"},
            "g2": {"contig": "chr", "start": 3000, "end": 3500, "strand": "+"},
        }

    def _design(self, genome, tid, monkeypatch):
        recs, coords = genome
        monkeypatch.setattr(design, "run_blast_batch",
                            lambda seqs, db, **k: {q: (1, [], True, None) for q in seqs})
        return design.process_single_target(tid, recs, coords, "db", False, "rc_tailed", False,
                                            num_candidates=10)

    def test_process_single_target(self, genome, monkeypatch):
        tid, pairs, err = self._design(genome, "g1", monkeypatch)
        assert err is None and pairs
        bed = pairs[0]["bed_row"].split("\t")
        assert " " not in bed[3] and "(" not in bed[3]
        assert int(bed[2]) - int(bed[1]) > 0

    def test_sort_key_is_numeric(self, genome, monkeypatch):
        _, pairs, _ = self._design(genome, "g1", monkeypatch)
        keys = [p["sort_key"] for p in pairs]
        assert all(isinstance(k[1], int) for k in keys)

    def test_packing_deterministic_with_seed(self, genome, monkeypatch):
        m = {}
        for t in ("g1", "g2"):
            _, pairs, _ = self._design(genome, t, monkeypatch)
            m[t] = sorted(pairs, key=lambda x: x["sort_key"])
        a = design.pack_targets_first_fit(m, False, 5, seed=3)
        b = design.pack_targets_first_fit(m, False, 5, seed=3)
        assert [[r["csv_row"]["target_id"] for r in p] for p in a] == \
               [[r["csv_row"]["target_id"] for r in p] for p in b]


class TestBiophysics:
    def test_invalid_sequence_is_not_reported_as_no_dimer(self):
        from core.biophysics import get_homodimer_dg
        assert get_homodimer_dg("ACGT" * 40) < -100  # primer3 rejects >60nt

    def test_heterodimer_symmetric(self):
        from core.biophysics import get_heterodimer_dg
        assert get_heterodimer_dg("GAGACCTACCAGCTCACCTC", "TCCACGGCAATCTCGATGAT") == \
               get_heterodimer_dg("TCCACGGCAATCTCGATGAT", "GAGACCTACCAGCTCACCTC")


class TestFeedback:
    def test_redesign_targets_mapped_to_gene_ids(self, tmp_path, monkeypatch):
        import argparse, json
        fb = {"target_recommendations": [
            {"target_id": "gyrB_amplicon_Strategy_0_2", "rel_mean": 0.0,
             "recommended_adjustment": 10.0, "action": "🔥 FLAG FOR REDESIGN"},
            {"target_id": "recA_amplicon_Strategy_0_0", "rel_mean": 1.0,
             "recommended_adjustment": 1.0, "action": "Maintain current concentration"}]}
        f = tmp_path / "fb.json"
        f.write_text(json.dumps(fb))
        seen = {}
        monkeypatch.setattr(design, "run_design_mode", lambda a: seen.update(
            targets=open(a.target_file).read().split(), rescue=a.rescue))
        args = argparse.Namespace(feedback=str(f), genome="g", gff="g",
                                  output_prefix=str(tmp_path / "o"), rescue=False)
        design.run_feedback_mode(args)
        assert seen == {"targets": ["gyrB"], "rescue": True}
        assert (tmp_path / "o_rebalancing_recipe.csv").exists()


class TestSimulatorAndFixture:
    def test_simulator_is_deterministic(self, tmp_path):
        import gzip
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(design.__file__)), "analysis"))
        import generate_simulated_fastq as sim
        genome = {"c": SeqRecord(Seq("ACGT" * 500), id="c")}
        iv = [("c", 100, 400, "t1")]
        outs = []
        for n in ("a", "b"):
            sim.generate_fastq_records(iv, genome, 20, 100, str(tmp_path / n), seed=5,
                                       non_uniform=True, thermal_bias=60.0)
            outs.append(gzip.open(str(tmp_path / n) + "_R1.fastq.gz", "rt").read())
        assert outs[0] == outs[1] and outs[0]

    def test_stable_offset_is_stable(self):
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(design.__file__)), "analysis"))
        import generate_simulated_fastq as sim
        assert sim.stable_offset("abc") == sim.stable_offset("abc")

    def test_fixture_builds(self, tmp_path):
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(design.__file__)), "scripts"))
        import make_test_fixture as fx
        out = fx.build("common_refs/ecoli_genome.cleaned.fna", str(tmp_path), coverage=5)
        for f in ("genome.fna", "targets.bed", "primers_trimming.fasta",
                  "FixtureSample_R1.fastq.gz", "FixtureSample_R2.fastq.gz"):
            assert os.path.exists(os.path.join(out, f))


class TestAnalysisUtils:
    def test_to_float_handles_picard_placeholders(self):
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(design.__file__)), "analysis"))
        from analysis_utils import to_float
        assert to_float("?") == 0.0 and to_float("") == 0.0 and to_float(None) == 0.0
        assert to_float("1.5") == 1.5

    def test_report_survives_undefined_metrics(self, tmp_path):
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(design.__file__)), "analysis"))
        import plot_coverage
        stats = [{'name': 'a', 'chrom': 'c', 'start': 1, 'end': 9, 'mean_coverage': 10.0}]
        plot_coverage.generate_html_report(
            {'MEAN_TARGET_COVERAGE': '10', 'PCT_AMPLIFIED_BASES': '1', 'FOLD_80_BASE_PENALTY': '?'},
            stats, str(tmp_path / "r.html"), str(tmp_path / "f.json"))
        assert (tmp_path / "f.json").exists()
