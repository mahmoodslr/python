from pathlib import Path
import csv
import matplotlib.pyplot as plt


DATA_FILE = Path(__file__).parent / "data" / "sequence.fasta"
NUCLEOTIDES = {"A", "T", "G", "C"}
START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}
BASE_DIR = Path(__file__).parent
DATA_FILE = BASE_DIR / "data" / "sequence.fasta"
REPORT_DIR = BASE_DIR / "reports"


# -------------------------
# Stage 1: Basic analysis
# -------------------------

def read_sequence(file_path=DATA_FILE):
    """Read a single DNA sequence from a FASTA file."""
    with open(file_path, "r", encoding="utf-8") as file:
        lines = file.readlines()
    title = None
    sequence_parts = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if title is not None:
                raise ValueError("Only one FASTA record is supported.")
            title = line[1:].strip()
        else:
            if title is None:
                raise ValueError("Invalid FASTA: missing header.")
            sequence_parts.append(line)
    if not title or not sequence_parts:
        raise ValueError("FASTA file is empty or incomplete.")
    sequence = clean_sequence("".join(sequence_parts))
    if not validate_sequence(sequence):
        raise ValueError("Invalid DNA sequence.")
    return title, sequence


def clean_sequence(sequence):
    """Normalize a DNA sequence."""
    return "".join(sequence.split()).upper()


def validate_sequence(sequence):
    """Check whether all characters are valid DNA nucleotides."""
    return bool(sequence) and set(sequence) <= NUCLEOTIDES


def calculate_length(sequence):
    return len(sequence)


def count_nucleotides(sequence):
    return {
        nucleotide: sequence.count(nucleotide)
        for nucleotide in "ATGC"
    }


def calculate_gc_content(sequence):
    if not sequence:
        return 0.0
    return (
        (sequence.count("G") + sequence.count("C"))
        / len(sequence)
        * 100
    )


def calculate_at_content(sequence):
    if not sequence:
        return 0.0
    return (
        (sequence.count("A") + sequence.count("T"))
        / len(sequence)
        * 100
    )


# -------------------------
# Stage 2: Biological analysis
# -------------------------

def get_complement(sequence):
    """Return the complementary DNA strand."""
    complement_table = str.maketrans("ATGC", "TACG")
    return sequence.translate(complement_table)


def get_reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    return get_complement(sequence)[::-1]


def get_reading_frames(sequence):
    """
    Return six reading frames:
    three forward frames and three reverse-strand frames.
    """
    reverse_sequence = get_reverse_complement(sequence)
    frames = {}
    for strand, dna in (
        ("+", sequence),
        ("-", reverse_sequence),
    ):
        for offset in range(3):
            frame_name = f"{strand}{offset + 1}"
            frames[frame_name] = [
                dna[i:i + 3]
                for i in range(offset, len(dna) - 2, 3)
            ]
    return frames


def find_codons(sequence):
    """
    Find start and stop codons in all six reading frames.

    Coordinates are 1-based and refer to the original
    DNA sequence, including for the reverse strand.
    """
    results = []
    length = len(sequence)
    for frame_name, codons in get_reading_frames(sequence).items():
        strand = frame_name[0]
        offset = int(frame_name[1]) - 1
        for index, codon in enumerate(codons):
            if codon not in STOP_CODONS and codon != START_CODON:
                continue
            oriented_start = offset + index * 3
            if strand == "+":
                start = oriented_start + 1
                end = oriented_start + 3
            else:
                start = length - oriented_start - 2
                end = length - oriented_start
            results.append({
                "codon": codon,
                "type": (
                    "start" if codon == START_CODON else "stop"
                ),
                "strand": strand,
                "frame": frame_name,
                "start": start,
                "end": end,
            })
    return results


def find_orfs(sequence):
    """
    Find candidate ORFs in all six reading frames.

    Each start codon is paired with the first downstream
    in-frame stop codon. Coordinates are 1-based and
    refer to the original sequence.
    """
    orfs = []
    length = len(sequence)
    for frame_name, codons in get_reading_frames(sequence).items():
        strand = frame_name[0]
        offset = int(frame_name[1]) - 1
        start_indices = []
        for index, codon in enumerate(codons):
            if codon == START_CODON:
                start_indices.append(index)
            elif codon in STOP_CODONS:
                for start_index in start_indices:
                    oriented_start = offset + start_index * 3
                    oriented_end = offset + index * 3 + 2
                    if strand == "+":
                        start = oriented_start + 1
                        end = oriented_end + 1
                    else:
                        start = length - oriented_end
                        end = length - oriented_start
                    orfs.append({
                        "strand": strand,
                        "frame": frame_name,
                        "start": start,
                        "end": end,
                        "length": oriented_end - oriented_start + 1,
                        "start_codon": START_CODON,
                        "stop_codon": codon,
                        "sequence": "".join(
                            codons[start_index:index + 1]
                        ),
                    })
                # A stop codon closes ORFs started earlier
                # in this same reading frame.
                start_indices = []
    return orfs


def find_pattern(sequence, pattern):
    """
    Find all occurrences of a DNA pattern.
    Overlapping matches are included.
    Return 1-based starting positions.
    """
    pattern = clean_sequence(pattern)
    if not pattern or not set(pattern) <= NUCLEOTIDES:
        raise ValueError("Pattern must contain only A, T, G and C.")
    positions = []
    start = 0
    while True:
        index = sequence.find(pattern, start)
        if index == -1:
            break
        positions.append(index + 1)
        start = index + 1
    return positions


def save_csv_report(
    title,
    sequence,
    codons,
    orfs,
    pattern="",
    pattern_positions=None,
):
    """Save a summary of the DNA analysis to a CSV file."""
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    counts = count_nucleotides(sequence)
    start_count = sum(
        1 for item in codons if item["type"] == "start"
    )
    stop_count = sum(
        1 for item in codons if item["type"] == "stop"
    )
    report = {
        "sequence_name": title,
        "sequence_length": calculate_length(sequence),
        "A_count": counts["A"],
        "T_count": counts["T"],
        "G_count": counts["G"],
        "C_count": counts["C"],
        "GC_content_percent": round(
            calculate_gc_content(sequence), 2
        ),
        "AT_content_percent": round(
            calculate_at_content(sequence), 2
        ),
        "start_codon_count": start_count,
        "stop_codon_count": stop_count,
        "candidate_orf_count": len(orfs),
        "searched_pattern": pattern,
        "pattern_match_count": len(pattern_positions or []),
        "pattern_positions": ";".join(
            str(position)
            for position in (pattern_positions or [])
        ),
    }
    report_file = REPORT_DIR / "analysis_report.csv"
    with open(
        report_file,
        "w",
        newline="",
        encoding="utf-8-sig",
    ) as file:
        writer = csv.DictWriter(
            file,
            fieldnames=report.keys(),
        )
        writer.writeheader()
        writer.writerow(report)
    return report_file


def plot_nucleotide_counts(sequence):
    """Create a bar chart of nucleotide counts."""
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    counts = count_nucleotides(sequence)
    nucleotides = list(counts.keys())
    values = list(counts.values())
    plt.figure(figsize=(8, 5))
    plt.bar(nucleotides, values)
    plt.title("DNA Nucleotide Counts")
    plt.xlabel("Nucleotide")
    plt.ylabel("Count")
    plt.tight_layout()
    output_file = REPORT_DIR / "nucleotide_counts.png"
    plt.savefig(output_file, dpi=150)
    plt.close()
    return output_file


def plot_gc_at_content(sequence):
    """Create a bar chart of GC and AT content."""
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    gc_content = calculate_gc_content(sequence)
    at_content = calculate_at_content(sequence)
    labels = ["GC Content", "AT Content"]
    values = [gc_content, at_content]
    plt.figure(figsize=(8, 5))
    plt.bar(labels, values)
    plt.title("DNA Base Composition")
    plt.ylabel("Percentage (%)")
    plt.ylim(0, 100)
    plt.tight_layout()
    output_file = REPORT_DIR / "gc_at_content.png"
    plt.savefig(output_file, dpi=150)
    plt.close()
    return output_file


# -------------------------
# Main program
# -------------------------

def main():
    try:
        title, sequence = read_sequence()
    except (OSError, ValueError) as error:
        print(f"Error: {error}")
        return
    print("\nDNA Sequence Analyzer")
    print("=" * 30)
    print(f"Sequence name: {title}")
    print(f"Sequence length: {calculate_length(sequence)}")
    print(f"Nucleotide counts: {count_nucleotides(sequence)}")
    print(f"GC content: {calculate_gc_content(sequence):.2f}%")
    print(f"AT content: {calculate_at_content(sequence):.2f}%")
    print("\nComplement:")
    print(get_complement(sequence))
    print("\nReverse complement:")
    print(get_reverse_complement(sequence))
    codons = find_codons(sequence)
    print("\nCodons:")
    for item in codons:
        print(item)
    print("\nReading frames:")
    for frame, frame_codons in get_reading_frames(sequence).items():
        print(f"{frame}: {' '.join(frame_codons)}")
    orfs = find_orfs(sequence)
    print("\nCandidate ORFs:")
    if orfs:
        for orf in orfs:
            print(orf)
    else:
        print("No candidate ORFs found.")
    pattern = input(
        "\nEnter a DNA pattern to search for "
        "(or press Enter to skip): "
    ).strip()
    pattern_positions = []
    if pattern:
        try:
            pattern_positions = find_pattern(sequence, pattern)
            print(f"Pattern positions (1-based): {pattern_positions}")
            print(f"Number of matches: {len(pattern_positions)}")
        except ValueError as error:
            print(f"Invalid pattern: {error}")
            return
    try:
        csv_file = save_csv_report(
            title,
            sequence,
            codons,
            orfs,
            pattern,
            pattern_positions,
        )
        counts_chart = plot_nucleotide_counts(sequence)
        content_chart = plot_gc_at_content(sequence)
    except OSError as error:
        print(f"Could not save reports: {error}")
        return
    print("\nReports generated successfully!")
    print(f"CSV report: {csv_file}")
    print(f"Nucleotide chart: {counts_chart}")
    print(f"GC/AT chart: {content_chart}")


if __name__ == "__main__":
    main()
