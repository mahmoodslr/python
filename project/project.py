import sys


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
    if len(sys.argv) != 3:
        sys.exit("Usage: python project.py <reference.fasta> <variant.fasta>")
    try:
        reference_id, reference_sequence = load_fasta(sys.argv[1])
        variant_id, variant_sequence = load_fasta(sys.argv[2])
        reference_sequence = validate_dna(reference_sequence)
        variant_sequence = validate_dna(variant_sequence)
        if len(reference_sequence) % 3 != 0:
            raise ValueError("Reference coding sequence length must be a multiple of 3.")
        aligned_reference, aligned_variant, score = needleman_wunsch(reference_sequence, variant_sequence)
        variants = detect_variants(aligned_reference, aligned_variant)
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
    print()
    print("Detected variants:")
    if not variants:
        print("No differences detected.")
    else:
        for number, variant in enumerate(variants, start=1):
            print(f"\nVariant {number}: {variant['type']}")
            for key, value in variant.items():
                if key != "type":
                    print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
