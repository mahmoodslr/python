import sys
import csv
from collections import Counter
from pathlib import Path


CODON_TABLE = {
    # Phenylalanine, Leucine
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",

    # Serine, Tyrosine, Stop, Cysteine, Tryptophan
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",

    # Leucine, Proline, Histidine, Glutamine, Arginine
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    # Isoleucine, Methionine, Threonine, Asparagine, Lysine
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",

    # Valine, Alanine, Aspartic acid, Glutamic acid, Glycine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


REPORT_FIELDS = [
    "reference_id",
    "variant_id",
    "reference_length_nt",
    "variant_length_nt",
    "net_length_change",
    "net_frame_phase_preserved",
    "reference_protein",
    "variant_protein",
    "protein_changed",
    "total_detected_variants",
    "variant_type_counts",
    "total_classified_effects",
    "effect_counts",
    "variant_type",
    "effect",
    "codon_number",
    "reference_codon",
    "variant_codon",
    "reference_amino_acid",
    "variant_amino_acid",
    "reference_positions",
    "reference_position",
    "variant_position",
    "reference_base",
    "variant_base",
    "after_reference_position",
    "variant_start",
    "variant_end",
    "inserted_sequence",
    "reference_start",
    "reference_end",
    "before_variant_position",
    "deleted_sequence",
    "length_nt",
    "notes",
    "reason"
]


def load_fasta(path):
    identifier= None
    sequence_parts= []
    with open(path, encoding="utf-8") as file:
        for line_number, raw_line in enumerate(file, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    raise ValueError("Expected exactly one FASTA record per file.")
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Missing sequence identifier on line {line_number}.")
                identifier = header.split()[0]
            else:
                if identifier is None:
                    raise ValueError(f"Sequence data before FASTA header , on line {line_number}.")
                sequence_parts.append(line)
    if identifier is None:
        raise ValueError("No FASTA header found.")
    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise ValueError("No DNA sequence found in FASTA record.")
    return identifier, sequence


def validate_dna(sequence):
    sequence = sequence.upper()
    if not sequence:
        raise ValueError("DNA sequence is empty.")
    invalid_symbols= sorted(set(sequence)-set("ACGT"))
    if invalid_symbols:
        symbols = ", ".join(repr(symbol) for symbol in invalid_symbols)
        raise ValueError(f"Invalid DNA symbol(s): {symbols}")
    return sequence


def translate_dna(sequence, *, allow_partial=False):
    sequence = validate_dna(sequence)
    remainder = len(sequence) % 3
    if remainder and not allow_partial:
        raise ValueError("DNA sequence length must be a multiple of 3.")
    complete_length = len(sequence) - remainder
    protein = []
    for i in range(0, complete_length, 3):
        codon = sequence[i:i + 3]
        amino_acid = CODON_TABLE[codon]
        protein.append(amino_acid)
        if amino_acid == "*":
            break
    return "".join(protein)


def classify_effects(reference_sequence, variant_sequence, variants):
    reference_sequence = validate_dna(reference_sequence)
    variant_sequence = validate_dna(variant_sequence)
    if len(reference_sequence) % 3 != 0:
        raise ValueError(
            "Reference coding sequence length must be a multiple of 3."
        )
    if not isinstance(variants, list):
        raise TypeError("variants must be a list.")
    known_types = {"substitution", "insertion", "deletion"}
    for variant in variants:
        if not isinstance(variant, dict):
            raise TypeError("Each variant must be a dictionary.")
        if variant.get("type") not in known_types:
            raise ValueError("Unknown variant type.")
    reference_protein = translate_dna(reference_sequence)
    variant_protein = translate_dna(
        variant_sequence,
        allow_partial=True
    )
    has_indels = any(
        variant["type"] in {"insertion", "deletion"}
        for variant in variants
    )
    effects = []
    notes = []
    # Case 1: substitutions without insertions or deletions.
    if not has_indels:
        if len(reference_sequence) != len(variant_sequence):
            raise ValueError(
                "Different sequence lengths require insertion "
                "or deletion events."
            )
        codons_affected = {}
        for variant in variants:
            if variant["type"] != "substitution":
                continue
            position = variant["reference_position"]
            codon_number = (position - 1) // 3 + 1
            codons_affected.setdefault(
                codon_number, []
            ).append(variant)
        for codon_number, codon_variants in sorted(
            codons_affected.items()
        ):
            start = (codon_number - 1) * 3
            end = start + 3
            reference_codon = reference_sequence[start:end]
            variant_codon = variant_sequence[start:end]
            if (
                len(reference_codon) != 3
                or len(variant_codon) != 3
            ):
                raise ValueError(
                    "Affected codon is incomplete."
                )
            reference_aa = CODON_TABLE[reference_codon]
            variant_aa = CODON_TABLE[variant_codon]
            if reference_aa == "*" and variant_aa == "*":
                effect = "stop_retained"
            elif reference_aa != "*" and variant_aa == "*":
                effect = "nonsense"
            elif reference_aa == "*" and variant_aa != "*":
                effect = "stop_loss"
            elif reference_aa == variant_aa:
                effect = "synonymous"
            else:
                effect = "missense"
            effects.append({
                "variant_type": "substitution",
                "effect": effect,
                "codon_number": codon_number,
                "reference_codon": reference_codon,
                "variant_codon": variant_codon,
                "reference_amino_acid": reference_aa,
                "variant_amino_acid": variant_aa,
                "reference_positions": [
                    item["reference_position"]
                    for item in codon_variants
                ]
            })
    # Case 2: insertions or deletions are present.
    else:
        for variant in variants:
            if variant["type"] == "insertion":
                inserted_sequence = variant["inserted_sequence"]
                length = len(inserted_sequence)
                effect = (
                    "frameshift_insertion"
                    if length % 3 != 0
                    else "in_frame_insertion"
                )
                effects.append({
                    "variant_type": "insertion",
                    "effect": effect,
                    "length_nt": length,
                    "inserted_sequence": inserted_sequence,
                    "after_reference_position":
                        variant["after_reference_position"]
                })
            elif variant["type"] == "deletion":
                deleted_sequence = variant["deleted_sequence"]
                length = len(deleted_sequence)
                effect = (
                    "frameshift_deletion"
                    if length % 3 != 0
                    else "in_frame_deletion"
                )
                effects.append({
                    "variant_type": "deletion",
                    "effect": effect,
                    "length_nt": length,
                    "deleted_sequence": deleted_sequence,
                    "reference_start": variant["reference_start"],
                    "reference_end": variant["reference_end"]
                })
            else:
                effects.append({
                    "variant_type": "substitution",
                    "effect": "requires_contextual_analysis",
                    "reference_position":
                        variant["reference_position"],
                    "reason": (
                        "Substitution consequences are not "
                        "classified independently when indels "
                        "are present."
                    )
                })
        notes.append(
            "Indels can alter the reading frame. Interpret "
            "protein-level differences using the full sequences."
        )
        notes.append(
            "The net length change alone does not reveal whether "
            "a temporary frameshift occurred between multiple indels."
        )
    if len(variant_sequence) % 3:
        notes.append(
            f"The variant ends with {len(variant_sequence) % 3} "
            "untranslated nucleotide(s) outside a complete codon."
        )
    return {
        "reference_protein": reference_protein,
        "variant_protein": variant_protein,
        "protein_changed": reference_protein != variant_protein,
        "reference_length_nt": len(reference_sequence),
        "variant_length_nt": len(variant_sequence),
        "net_length_change": (
            len(variant_sequence) - len(reference_sequence)
        ),
        "net_frame_phase_preserved": (
            (len(variant_sequence) - len(reference_sequence)) % 3 == 0
        ),
        "effects": effects,
        "notes": notes
    }


def summarize_results(variants, effect_report):
    variant_type_counts = dict(sorted(Counter(item["type"] for item in variants).items()))
    effects = effect_report["effects"]
    effect_counts = dict(sorted(Counter(item["effect"] for item in effects).items()))
    return {
        "total_detected_variants": len(variants),
        "variant_type_counts": variant_type_counts,
        "total_classified_effects": len(effects),
        "effect_counts": effect_counts,
        "protein_changed": effect_report["protein_changed"],
        "net_length_change": effect_report["net_length_change"],
        "net_frame_phase_preserved":
            effect_report["net_frame_phase_preserved"]
    }


def write_report_csv(
    output_path,
    reference_id,
    variant_id,
    variants,
    effect_report,
    summary
):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    effects = effect_report["effects"]
    if not effects:
        effects = [{
            "variant_type": "none",
            "effect": "no_change"
        }]
    variant_counts_text = "; ".join(
        f"{key}: {value}"
        for key, value in summary["variant_type_counts"].items()
    )
    effect_counts_text = "; ".join(
        f"{key}: {value}"
        for key, value in summary["effect_counts"].items()
    )
    with output_path.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as file:
        writer = csv.DictWriter(
            file,
            fieldnames=REPORT_FIELDS
        )
        writer.writeheader()
        for effect in effects:
            row = {
                "reference_id": reference_id,
                "variant_id": variant_id,
                "reference_length_nt":
                    effect_report["reference_length_nt"],
                "variant_length_nt":
                    effect_report["variant_length_nt"],
                "net_length_change":
                    summary["net_length_change"],
                "net_frame_phase_preserved":
                    summary["net_frame_phase_preserved"],
                "reference_protein":
                    effect_report["reference_protein"],
                "variant_protein":
                    effect_report["variant_protein"],
                "protein_changed":
                    effect_report["protein_changed"],
                "total_detected_variants":
                    summary["total_detected_variants"],
                "variant_type_counts": variant_counts_text,
                "total_classified_effects":
                    summary["total_classified_effects"],
                "effect_counts": effect_counts_text,
                "notes": " | ".join(effect_report["notes"])
            }
            row.update(effect)
            if isinstance(row.get("reference_positions"), list):
                row["reference_positions"] = "; ".join(
                    map(str, row["reference_positions"])
                )
            writer.writerow(row)
    return output_path


def needleman_wunsch(
    reference,
    variant,
    match_score=1,
    mismatch_score=-1,
    gap_score=-1
):
    if not isinstance(reference, str) or not isinstance(variant, str):
        raise TypeError("Sequences must be strings.")
    reference = reference.upper()
    variant = variant.upper()
    invalid_reference = set(reference) - set("ACGT")
    invalid_variant = set(variant) - set("ACGT")
    if invalid_reference or invalid_variant:
        raise ValueError("Sequences must contain only A, C, G, and T.")
    n = len(reference)
    m = len(variant)
    # 1. Create the scoring matrix.
    matrix = [[0]* (m + 1) for _ in range(n + 1)]
    # 2. Initialize the first column and row.
    for i in range(1, n + 1):
        matrix[i][0] = matrix[i - 1][0] + gap_score
    for j in range(1, m + 1):
        matrix[0][j] = matrix[0][j - 1] + gap_score
    # 3. Calculate the best score for every cell.
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if reference[i - 1] == variant[j - 1]:
                substitution_score = match_score
            else:
                substitution_score = mismatch_score
            diagonal = matrix[i - 1][j - 1] + substitution_score
            up = matrix[i - 1][j] + gap_score
            left = matrix[i][j - 1] + gap_score
            matrix[i][j] = max(diagonal, up, left)
    # 4. Reconstruct the alignment.
    aligned_reference = []
    aligned_variant = []
    i = n
    j = m
    while i>0 or j>0:
        if i>0 and j>0:
            if reference[i - 1] == variant[j - 1]:
                substitution_score = match_score
            else:
                substitution_score = mismatch_score
            if matrix[i][j] == (matrix[i - 1][j - 1] + substitution_score):
                aligned_reference.append(reference[i - 1])
                aligned_variant.append(variant[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and matrix[i][j] == (
            matrix[i - 1][j] + gap_score
        ):
            aligned_reference.append(reference[i - 1])
            aligned_variant.append("-")
            i -= 1
            continue
        if j > 0 and matrix[i][j] == (
            matrix[i][j - 1] + gap_score
        ):
            aligned_reference.append("-")
            aligned_variant.append(variant[j - 1])
            j -= 1
            continue
        raise RuntimeError("Could not reconstruct the alignment.")
    aligned_reference.reverse()
    aligned_variant.reverse()
    return (
        "".join(aligned_reference),
        "".join(aligned_variant),
        matrix[n][m]
    )


def detect_variants(aligned_reference, aligned_variant):
    if not isinstance(aligned_reference, str) or not isinstance(aligned_variant, str):
        raise TypeError("Aligned sequences must be strings.")
    aligned_reference = aligned_reference.upper()
    aligned_variant = aligned_variant.upper()
    if len(aligned_reference) != len(aligned_variant):
        raise ValueError("Aligned sequences must have equal lengths.")
    allowed_symbols = set("ACGT-")
    invalid_symbols = (set(aligned_reference) | set(aligned_variant)) - allowed_symbols
    if invalid_symbols:
        symbols = ", ".join(sorted(invalid_symbols))
        raise ValueError(f"Invalid alignment symbols: {symbols}")
    if any(ref == "-" and var == "-" for ref, var in zip(aligned_reference, aligned_variant)):
        raise ValueError("An alignment column cannot contain two gaps.")
    variants = []
    reference_position = 0
    variant_position = 0
    i = 0
    while i < len(aligned_reference):
        ref = aligned_reference[i]
        var = aligned_variant[i]
        # Case 1: Insertion relative to the reference.
        if ref == "-":
            insertion_after_reference = reference_position
            variant_start = variant_position + 1
            inserted_bases = []
            while (
                i < len(aligned_reference)
                and aligned_reference[i] == "-"
                and aligned_variant[i] != "-"
            ):
                inserted_bases.append(aligned_variant[i])
                variant_position += 1
                i += 1
            variants.append({
                "type": "insertion",
                "after_reference_position": insertion_after_reference,
                "variant_start": variant_start,
                "variant_end": variant_position,
                "inserted_sequence": "".join(inserted_bases)
            })
        # Case 2: Deletion relative to the reference.
        elif var == "-":
            reference_start = reference_position + 1
            variant_boundary = variant_position
            deleted_bases = []
            while (
                i < len(aligned_reference)
                and aligned_reference[i] != "-"
                and aligned_variant[i] == "-"
            ):
                deleted_bases.append(aligned_reference[i])
                reference_position += 1
                i += 1

            variants.append({
                "type": "deletion",
                "reference_start": reference_start,
                "reference_end": reference_position,
                "before_variant_position": variant_boundary,
                "deleted_sequence": "".join(deleted_bases)
            })
        # Case 3: Both sequences have a nucleotide.
        else:
            reference_position += 1
            variant_position += 1
            if ref != var:
                variants.append({
                    "type": "substitution",
                    "reference_position": reference_position,
                    "variant_position": variant_position,
                    "reference_base": ref,
                    "variant_base": var
                })
            i += 1
    return variants


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit("Usage: python project.py <reference.fasta> <variant.fasta> [report.csv]")
    output_path = sys.argv[3] if len(sys.argv) == 4 else "analysis.csv"
    try:
        input_paths = {
            Path(sys.argv[1]).resolve(),
            Path(sys.argv[2]).resolve()
        }
        if Path(output_path).resolve() in input_paths:
            raise ValueError("Output report must not overwrite an input FASTA file.")
        reference_id, reference_sequence = load_fasta(sys.argv[1])
        variant_id, variant_sequence = load_fasta(sys.argv[2])
        reference_sequence = validate_dna(reference_sequence)
        variant_sequence = validate_dna(variant_sequence)
        if len(reference_sequence) % 3 != 0:
            raise ValueError("Reference coding sequence length must be a multiple of 3.")
        aligned_reference, aligned_variant, score = needleman_wunsch(reference_sequence, variant_sequence)
        variants = detect_variants(aligned_reference, aligned_variant)
        effect_report = classify_effects(reference_sequence, variant_sequence, variants)
        summary = summarize_results(variants, effect_report) 
        report_path = write_report_csv(
            output_path,
            reference_id,
            variant_id,
            variants,
            effect_report,
            summary
        )
        # reference_protein = translate_dna(reference_sequence)
        # variant_protein = translate_dna(variant_sequence, allow_partial=True)
    except FileNotFoundError as error:
        sys.exit(f"File not found: {error.filename}")
    except (OSError, ValueError) as error:
        sys.exit(f"Input error: {error}")
    print("MutationScope - Sequence Analysis")
    print(f"Reference: {reference_id} ({len(reference_sequence)} nt)")
    print(f"Variant: {variant_id} ({len(variant_sequence)} nt)")
    print("\nAlignment:")
    print(f"Reference: {aligned_reference}")
    print(f"Variant:   {aligned_variant}")
    print(f"Alignment score: {score}")
    print("\nDetected variants:")
    if not variants:
        print("No differences detected.")
    else:
        for number, variant in enumerate(variants, start=1):
            print(f"\nVariant {number}: {variant['type']}")
            for key, value in variant.items():
                if key != "type":
                    print(f"  {key}: {value}")
    print("\nProtein translation:")
    # print(f"Reference protein: {reference_protein}")
    # print(f"Variant protein:   {variant_protein}")
    print("\nProtein-level analysis:")
    print(f"Reference protein: {effect_report['reference_protein']}")
    print(f"Variant protein:   {effect_report['variant_protein']}")
    print(f"Protein changed:   {effect_report['protein_changed']}")
    print("\nEffect classification:")
    if not effect_report["effects"]:
        print("No coding effects identified.")
    else:
        for number, effect in enumerate(
            effect_report["effects"], start=1
        ):
            print(f"\nEffect {number}: {effect['effect']}")
            for key, value in effect.items():
                if key != "effect":
                    print(f"  {key}: {value}")
    if effect_report["notes"]:
        print("\nNotes:")
        for note in effect_report["notes"]:
            print(f"- {note}")
    print(f"\nReport saved to: {report_path}")


if __name__ == "__main__":
    main()
