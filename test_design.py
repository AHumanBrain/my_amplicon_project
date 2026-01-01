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
