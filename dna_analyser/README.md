# DNA Sequence Analyzer

A command-line Python program that reads one DNA sequence from a FASTA file,
validates it, analyzes it, and saves a CSV report and two charts.

## Features

**Stage 1: Basic analysis**

- Reads the title and sequence from a FASTA file
- Cleans the sequence (removes spaces and line breaks, converts to uppercase)
- Validates the sequence (only `A`, `T`, `G` and `C` are accepted)
- Calculates sequence length, nucleotide counts, and GC / AT percentages
- Handles errors: missing file, empty file, missing header, several records, invalid characters

**Stage 2: Biological analysis**

- Complement and reverse complement
- Six reading frames (`+1`, `+2`, `+3`, `-1`, `-2`, `-3`)
- Start codon (`ATG`) and stop codons (`TAA`, `TAG`, `TGA`) with positions
- Candidate ORFs with strand, frame, start, end, length and sequence
- Search for a custom DNA pattern (overlapping matches are included)

**Stage 3: Reports and charts**

- Summary printed in the terminal
- CSV report saved to `reports/analysis_report.csv`
- Bar chart of `A`, `T`, `G`, `C` counts saved to `reports/nucleotide_counts.png`
- Bar chart of GC and AT percentages saved to `reports/gc_at_content.png`

## Project structure

```
dna_analyzer/
├── project.py          # main program and all analysis functions
├── test_project.py     # automated tests
├── README.md
├── data/
│   └── sequence.fasta  # input file
└── reports/
    ├── analysis_report.csv
    ├── nucleotide_counts.png
    └── gc_at_content.png
```

The `reports/` folder is created automatically when the program runs.

## Requirements

- Python 3.8 or newer
- `matplotlib` (for the charts)
- `pytest` (only for running the tests)

```
pip install matplotlib pytest
```

## Input format

Put one FASTA record in `data/sequence.fasta`:

```
>my_sequence
ATGAAATAGCCATGCGTACGTTAGCTAGC
ATGGGCTAACGATCGATCG
```

The sequence may be split over several lines and may use lowercase letters.

## Usage

```
python project.py
```

The program prints the analysis, then asks:

```
Enter a DNA pattern to search for (or press Enter to skip):
```

Type a pattern such as `ATG` to find all its positions, or press Enter to skip.
Then the CSV report and the charts are saved in `reports/`.

## CSV report columns

`sequence_name`, `sequence_length`, `A_count`, `T_count`, `G_count`, `C_count`,
`GC_content_percent`, `AT_content_percent`, `start_codon_count`,
`stop_codon_count`, `candidate_orf_count`, `searched_pattern`,
`pattern_match_count`, `pattern_positions`

Multiple pattern positions are separated by `;`.

## Running the tests

```
pytest
```

The tests cover the analysis functions, FASTA error handling, ORF coordinates
on both strands, pattern search, and the CSV and chart files. Report and chart
tests use a temporary folder, so your real `reports/` folder is not changed.

## Design notes and limitations

- All positions are **1-based** and refer to the **original** sequence,
  also for the reverse strand.
- An ORF is a start codon (`ATG`) paired with the first in-frame stop codon
  after it. Sequences with a start codon but no stop codon are not reported.
- If several `ATG` codons appear before the same stop codon, each one is
  reported as its own (nested) ORF.
- **A candidate ORF is not proof of a real gene.** In a random sequence, short
  ORFs appear by chance. Real gene prediction needs more evidence.
- Only one FASTA record is supported, and only `A`, `T`, `G`, `C`
  (no `N`, `U` or IUPAC ambiguity codes).
- Pattern search looks at the forward strand only.
- If the pattern entered by the user is invalid, the program prints an error
  and stops before saving the reports.
