import pytest
from project import load_fasta, validate_dna, needleman_wunsch, detect_variants, translate_dna, classify_effects,summarize_results,write_report_csv, CODON_TABLE
import csv


def test_two_substitutions_are_not_misclassified_as_indels():
    aligned_reference, aligned_variant, score = needleman_wunsch(
        "ACGT",
        "AGCT"
    )
    assert "-" not in aligned_reference
    assert "-" not in aligned_variant
    variants = detect_variants(aligned_reference, aligned_variant)
    assert [variant["type"] for variant in variants] == [
        "substitution",
        "substitution"
    ]
    assert [
        (
            variant["reference_position"],
            variant["reference_base"],
            variant["variant_base"]
        )
        for variant in variants
    ] == [
        (2, "C", "G"),
        (3, "G", "C")
    ]


def test_single_nucleotide_insertion_is_still_detected():
    aligned_reference, aligned_variant, score = needleman_wunsch(
        "ACGT",
        "ACAGT"
    )
    variants = detect_variants(aligned_reference, aligned_variant)
    assert len(variants) == 1
    assert variants[0]["type"] == "insertion"
    assert variants[0]["inserted_sequence"] == "A"


def test_single_nucleotide_deletion_is_still_detected():
    aligned_reference, aligned_variant, score = needleman_wunsch(
        "ACAGT",
        "ACGT"
    )
    variants = detect_variants(aligned_reference, aligned_variant)
    assert len(variants) == 1
    assert variants[0]["type"] == "deletion"
    assert variants[0]["deleted_sequence"] == "A"


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


def test_detect_variants_no_change():
    assert detect_variants("ACGT", "ACGT") == []


def test_detect_variants_substitution():
    result = detect_variants("ATG", "ATA")
    assert result == [{
        "type": "substitution",
        "reference_position": 3,
        "variant_position": 3,
        "reference_base": "G",
        "variant_base": "A"
    }]


def test_detect_variants_insertion():
    result = detect_variants("AC-GT", "ACTGT")
    assert result == [{
        "type": "insertion",
        "after_reference_position": 2,
        "variant_start": 3,
        "variant_end": 3,
        "inserted_sequence": "T"
    }]


def test_detect_variants_deletion():
    result = detect_variants("ACTTGT", "AC--GT")
    assert result == [{
        "type": "deletion",
        "reference_start": 3,
        "reference_end": 4,
        "before_variant_position": 2,
        "deleted_sequence": "TT"
    }]


def test_detect_variants_rejects_unequal_alignment_lengths():
    with pytest.raises(ValueError, match="equal lengths"):
        detect_variants("ACGT", "ACG")


def test_detect_variants_rejects_invalid_symbols():
    with pytest.raises(ValueError, match="Invalid alignment symbols"):
        detect_variants("ACXT", "ACGT")


def test_detect_variants_rejects_double_gap():
    with pytest.raises(ValueError, match="two gaps"):
        detect_variants("A-C", "A-C")


def test_reference_variant_workflow():
    reference = "ATGGAACCTTGG"
    variant = "ATGGTACCTTGG"
    aligned_reference, aligned_variant, score = needleman_wunsch(
        reference, variant
    )
    variants = detect_variants(aligned_reference, aligned_variant)
    assert aligned_reference == reference
    assert aligned_variant == variant
    assert score == 10
    assert variants == [
        {
            "type": "substitution",
            "reference_position": 5,
            "variant_position": 5,
            "reference_base": "A",
            "variant_base": "T"
        }
    ]


def test_codon_table_has_64_codons():
    assert len(CODON_TABLE) == 64


def test_translate_reference_sequence():
    assert translate_dna("ATGGAACCTTGG") == "MEPW"


def test_translate_variant_sequence():
    assert translate_dna("ATGGTACCTTGG") == "MVPW"


def test_translate_stops_at_first_stop():
    assert translate_dna("ATGTAAACG") == "M*"


def test_translate_rejects_incomplete_sequence():
    with pytest.raises(ValueError, match="multiple of 3"):
        translate_dna("ATGG")


def test_translate_allows_incomplete_terminal_codon():
    assert translate_dna("ATGG", allow_partial=True) == "M"


def analyze_pair(reference, variant):
    aligned_reference, aligned_variant, _ = needleman_wunsch(
        reference, variant
    )
    variants = detect_variants(
        aligned_reference, aligned_variant
    )
    return classify_effects(reference, variant, variants)


def test_classify_synonymous():
    result = analyze_pair("ATGGAACCTTGG", "ATGGAGCCTTGG")
    assert result["reference_protein"] == "MEPW"
    assert result["variant_protein"] == "MEPW"
    assert result["protein_changed"] is False
    assert result["effects"][0]["effect"] == "synonymous"


def test_classify_missense():
    result = analyze_pair("ATGGAACCTTGG", "ATGGTACCTTGG")
    assert result["reference_protein"] == "MEPW"
    assert result["variant_protein"] == "MVPW"
    assert result["protein_changed"] is True
    assert result["effects"][0]["effect"] == "missense"


def test_classify_nonsense():
    result = analyze_pair("ATGGAACCTTGG", "ATGTAACCTTGG")
    assert result["reference_protein"] == "MEPW"
    assert result["variant_protein"] == "M*"
    assert result["effects"][0]["effect"] == "nonsense"


def test_classify_stop_loss():
    result = analyze_pair("ATGTAA", "ATGCAA")
    assert result["reference_protein"] == "M*"
    assert result["variant_protein"] == "MQ"
    assert result["effects"][0]["effect"] == "stop_loss"


def test_classify_stop_retained():
    result = analyze_pair("ATGTAA", "ATGTAG")
    assert result["reference_protein"] == "M*"
    assert result["variant_protein"] == "M*"
    assert result["effects"][0]["effect"] == "stop_retained"


def test_classify_frameshift_insertion():
    reference = "ATGGAACCTTGG"
    variant = "ATGAGAACCTTGG"
    aligned_reference = "ATG-GAACCTTGG"
    aligned_variant = "ATGAGAACCTTGG"
    variants = detect_variants(
        aligned_reference, aligned_variant
    )
    result = classify_effects(reference, variant, variants)
    assert result["effects"][0]["effect"] == "frameshift_insertion"
    assert result["effects"][0]["length_nt"] == 1
    assert result["protein_changed"] is True


def test_classify_in_frame_deletion():
    reference = "ATGCCCGAA"
    variant = "ATGGAA"
    aligned_reference = "ATGCCCGAA"
    aligned_variant = "ATG---GAA"
    variants = detect_variants(
        aligned_reference, aligned_variant
    )
    result = classify_effects(reference, variant, variants)
    assert result["effects"][0]["effect"] == "in_frame_deletion"
    assert result["effects"][0]["length_nt"] == 3
    assert result["reference_protein"] == "MPE"
    assert result["variant_protein"] == "ME"


def test_summarize_results():
    reference = "ATGGAACCTTGG"
    variant = "ATGGTACCTTGG"
    aligned_reference, aligned_variant, _ = needleman_wunsch(
        reference, variant
    )
    variants = detect_variants(aligned_reference, aligned_variant)
    effect_report = classify_effects(
        reference, variant, variants
    )
    summary = summarize_results(variants, effect_report)
    assert summary["total_detected_variants"] == 1
    assert summary["variant_type_counts"] == {
        "substitution": 1
    }
    assert summary["effect_counts"] == {
        "missense": 1
    }
    assert summary["protein_changed"] is True


def test_write_report_csv(tmp_path):
    reference = "ATGGAACCTTGG"
    variant = "ATGGTACCTTGG"
    aligned_reference, aligned_variant, _ = needleman_wunsch(
        reference, variant
    )
    variants = detect_variants(aligned_reference, aligned_variant)
    effect_report = classify_effects(
        reference, variant, variants
    )
    summary = summarize_results(variants, effect_report)
    output_path = tmp_path / "analysis.csv"
    write_report_csv(
        output_path,
        "reference",
        "variant",
        variants,
        effect_report,
        summary
    )
    with output_path.open(
        newline="",
        encoding="utf-8"
    ) as file:
        rows = list(csv.DictReader(file))
    assert len(rows) == 1
    assert rows[0]["reference_id"] == "reference"
    assert rows[0]["variant_id"] == "variant"
    assert rows[0]["effect"] == "missense"
    assert rows[0]["reference_positions"] == "5"
    assert rows[0]["reference_protein"] == "MEPW"
    assert rows[0]["variant_protein"] == "MVPW"


def test_substitution_before_insertion_is_classified():
    reference = "ATGGAACCTTGG"
    variant = "ATGGTACCTTGGA"
    variants = [
        {
            "type": "substitution",
            "reference_position": 5,
            "variant_position": 5,
            "reference_base": "A",
            "variant_base": "T",
        },
        {
            "type": "insertion",
            "after_reference_position": 12,
            "variant_start": 13,
            "variant_end": 13,
            "inserted_sequence": "A",
        },
    ]
    report = classify_effects(reference, variant, variants)
    substitution = next(
        effect
        for effect in report["effects"]
        if effect["variant_type"] == "substitution"
    )
    assert substitution["effect"] == "missense"


def test_substitution_during_frameshift_requires_context():
    reference = "ATGGAACCTTGG"
    variant = "ATGAGTACCTTGG"
    variants = [
        {
            "type": "insertion",
            "after_reference_position": 3,
            "variant_start": 4,
            "variant_end": 4,
            "inserted_sequence": "A",
        },
        {
            "type": "substitution",
            "reference_position": 5,
            "variant_position": 6,
            "reference_base": "A",
            "variant_base": "T",
        },
    ]
    report = classify_effects(reference, variant, variants)
    substitution = next(
        effect
        for effect in report["effects"]
        if effect["variant_type"] == "substitution"
    )
    assert substitution["effect"] == "requires_contextual_analysis"


def test_compensating_indels_restore_frame_phase():
    reference = "ATGGAACCTTGG"
    variant = "ATGAGAACTTGG"
    variants = [
        {
            "type": "insertion",
            "after_reference_position": 3,
            "variant_start": 4,
            "variant_end": 4,
            "inserted_sequence": "A",
        },
        {
            "type": "deletion",
            "reference_start": 7,
            "reference_end": 7,
            "before_variant_position": 6,
            "deleted_sequence": "C",
        },
    ]
    report = classify_effects(reference, variant, variants)
    insertion = next(
        effect
        for effect in report["effects"]
        if effect["variant_type"] == "insertion"
    )
    deletion = next(
        effect
        for effect in report["effects"]
        if effect["variant_type"] == "deletion"
    )
    assert insertion["frame_phase_after_event"] == 1
    assert deletion["frame_phase_after_event"] == 0
    assert deletion["frame_restored_after_event"] is True
    assert report["net_frame_phase_preserved"] is True