import csv
import matplotlib
matplotlib.use("Agg")
import pytest
import project


def test_clean_sequence():
    assert project.clean_sequence(" atg c\n") == "ATGC"


def test_validate_sequence():
    assert project.validate_sequence("ATGC") is True
    assert project.validate_sequence("ATGN") is False
    assert project.validate_sequence("") is False


def test_calculate_length():
    assert project.calculate_length("ATGCA") == 5


def test_count_nucleotides():
    result = project.count_nucleotides("AATGCC")
    assert result == {"A": 2, "T": 1, "G": 1, "C": 2}


def test_gc_content():
    assert project.calculate_gc_content("ATGC") == 50
    assert project.calculate_gc_content("GGCC") == 100
    assert project.calculate_gc_content("") == 0


def test_at_content():
    assert project.calculate_at_content("ATGC") == 50
    assert project.calculate_at_content("AATT") == 100
    assert project.calculate_at_content("") == 0


def test_read_sequence_valid(tmp_path):
    file = tmp_path / "test.fasta"
    file.write_text(">seq1\nATGC\natgc\n")
    title, sequence = project.read_sequence(file)
    assert title == "seq1"
    assert sequence == "ATGCATGC"


def test_read_sequence_missing_file(tmp_path):
    with pytest.raises(OSError):
        project.read_sequence(tmp_path / "no_file.fasta")


def test_read_sequence_empty_file(tmp_path):
    file = tmp_path / "empty.fasta"
    file.write_text("")
    with pytest.raises(ValueError):
        project.read_sequence(file)


def test_read_sequence_invalid_dna(tmp_path):
    file = tmp_path / "bad.fasta"
    file.write_text(">bad\nATGXYZ\n")
    with pytest.raises(ValueError):
        project.read_sequence(file)


def test_complement():
    assert project.get_complement("ATGC") == "TACG"


def test_reverse_complement():
    assert project.get_reverse_complement("ATGC") == "GCAT"


def test_reading_frames():
    frames = project.get_reading_frames("ATGAAATAG")
    assert len(frames) == 6
    assert frames["+1"] == ["ATG", "AAA", "TAG"]


def test_find_codons():
    codons = project.find_codons("ATGAAATAG")
    types = [item["type"] for item in codons]
    assert "start" in types
    assert "stop" in types


def test_find_orfs_simple():
    orfs = project.find_orfs("ATGAAATAG")
    assert len(orfs) == 1
    assert orfs[0]["frame"] == "+1"
    assert orfs[0]["start"] == 1
    assert orfs[0]["end"] == 9
    assert orfs[0]["length"] == 9
    assert orfs[0]["stop_codon"] == "TAG"


def test_find_orfs_reverse_strand():
    orfs = project.find_orfs("CTATTTCAT")
    assert len(orfs) == 1
    assert orfs[0]["strand"] == "-"
    assert orfs[0]["start"] == 1
    assert orfs[0]["end"] == 9


def test_find_orfs_no_stop_codon():
    assert project.find_orfs("ATGAAACCC") == []


def test_find_orfs_no_start_codon():
    assert project.find_orfs("CCCAAATAG") == []


def test_find_pattern():
    assert project.find_pattern("ATGCATGC", "ATG") == [1, 5]


def test_find_pattern_overlapping():
    assert project.find_pattern("ATATAT", "ATA") == [1, 3]


def test_find_pattern_not_found():
    assert project.find_pattern("ATGC", "GGG") == []


def test_find_pattern_invalid():
    with pytest.raises(ValueError):
        project.find_pattern("ATGC", "XYZ")


def test_save_csv_report(tmp_path, monkeypatch):
    monkeypatch.setattr(project, "REPORT_DIR", tmp_path)
    sequence = "ATGAAATAG"
    codons = project.find_codons(sequence)
    orfs = project.find_orfs(sequence)
    report_file = project.save_csv_report("demo", sequence, codons, orfs)
    assert report_file.exists()
    with open(report_file, newline="", encoding="utf-8-sig") as file:
        row = next(csv.DictReader(file))
    assert row["sequence_name"] == "demo"
    assert row["sequence_length"] == "9"
    assert row["candidate_orf_count"] == "1"


def test_plot_nucleotide_counts(tmp_path, monkeypatch):
    monkeypatch.setattr(project, "REPORT_DIR", tmp_path)
    output_file = project.plot_nucleotide_counts("ATGCATGC")
    assert output_file.exists()


def test_plot_gc_at_content(tmp_path, monkeypatch):
    monkeypatch.setattr(project, "REPORT_DIR", tmp_path)
    output_file = project.plot_gc_at_content("ATGCATGC")
    assert output_file.exists()