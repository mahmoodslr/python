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



def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python project.py, <reference.fasta> <variant.fasta>")
    try:
        reference_id, reference_sequence = load_fasta(sys.argv[1])
        variant_id, variant_sequence = load_fasta(sys.argv[2])
        reference_sequence = validate_dna(reference_sequence)
        variant_sequence = validate_dna(variant_sequence)
        if len(reference_sequence) % 3 != 0:
            raise ValueError("Reference coding sequence length must be a multiple of 3.")
    except FileNotFoundError as error:
        sys.exit(f"File not found: {error.filename}")
    except (OSError, ValueError) as error:
        sys.exit(f"Input error: {error}")

    print("MutationScope - Input Validation")
    print(f"Reference: {reference_id} ({len(reference_sequence)} nt)")
    print(f"Variant: {variant_id} ({len(variant_sequence)} nt)")
    print("Both DNA sequences passed basic validation.")


if __name__ == "__main__":
    main()
