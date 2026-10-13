import pytest
from project import load_fasta, validate_dna,needleman_wunsch


def test_validate_dna_accepts_canonical_bases():
    assert validate_dna("atgc") == "ATGC"


def test_validate_dna_rejects_empty_sequence():
    with pytest.raises(ValueError, match="empty"):
        validate_dna("")


def test_validate_dna_rejects_invalid_symbols():
    with pytest.raises(ValueError, match="Invalid DNA symbol"):
        validate_dna("ATGX")


def test_load_fasta_joins_wrapped_lines(tmp_path):
    path = tmp_path / "reference.fasta"
    path.write_text(
        ">reference coding sequence\nATGGAA\nCCTTGG\n",
        encoding="utf-8"
    )
    identifier, sequence = load_fasta(path)
    assert identifier == "reference"
    assert sequence == "ATGGAACCTTGG"


def test_load_fasta_rejects_sequence_before_header(tmp_path):
    path = tmp_path / "invalid.fasta"
    path.write_text("ATGGAACCTTGG\n", encoding="utf-8")
    with pytest.raises(ValueError, match="before FASTA header"):
        load_fasta(path)


def test_load_fasta_rejects_multiple_records(tmp_path):
    path = tmp_path / "multiple.fasta"
    path.write_text(
        ">reference\nATG\n>variant\nGAA\n",
        encoding="utf-8"
    )

    with pytest.raises(ValueError, match="exactly one FASTA record"):
        load_fasta(path)



def test_identical_sequences():
    ref, var, score = needleman_wunsch("ACGT", "ACGT")

    assert ref == "ACGT"
    assert var == "ACGT"
    assert score == 4


def test_insertion():
    ref, var, score = needleman_wunsch("ACGT", "ACTGT")

    assert ref == "AC-GT"
    assert var == "ACTGT"
    assert score == 3


def test_deletion():
    ref, var, score = needleman_wunsch("ACTGT", "ACGT")

    assert ref == "ACTGT"
    assert var == "AC-GT"
    assert score == 3


def test_substitution():
    ref, var, score = needleman_wunsch("ATG", "ATA")

    assert ref == "ATG"
    assert var == "ATA"
    assert score == 1
