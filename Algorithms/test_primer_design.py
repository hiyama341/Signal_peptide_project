import pytest
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from primer_design import generate_primers

def test_generate_primers():
    # Define test inputs and expected outputs
    nucleotide_sequences = ["ATGGCAGTTCGT", "AGTACGTGCTGA", SeqRecord(Seq("GATCTGCATGAC"), id="test_seq")]
    up_homology_arm = "ACGT"
    down_homology_arm = "TTGAC"
    nucleotide_homology = 3

    expected_output = pd.DataFrame({
        "Forward primer": ["actgccATACGT", "cacgtaCTACGT", "tgcagaTCACGT"],
        "Forward primer length": [12, 12, 12],
        "Reverse primer": ["ggcagtTCGTTTGAC", "tacgtgCTGATTGAC", "tctgcaTGACTTGAC"],
        "Reverse primer length": [15, 15, 15]
    })

    # Ensure expected output is generated
    output = generate_primers(nucleotide_sequences, up_homology_arm, down_homology_arm, nucleotide_homology)
    assert output.equals(expected_output)

    # Ensure TypeError is raised if nucleotide_sequences contains an invalid type
    with pytest.raises(TypeError):
        generate_primers([1, 2, 3], up_homology_arm, down_homology_arm, nucleotide_homology)

    # Ensure ValueError is raised if nucleotide_sequences contains non-nucleotide characters
    with pytest.raises(ValueError, match="Sequence .* contains non-nucleotide characters"):
        generate_primers(["ATGCN", SeqRecord(Seq("ATGCN"), id="test_seq")], up_homology_arm, down_homology_arm, nucleotide_homology)

     # Stress test with a large number of input sequences
    with open("../Data/FungiDB-62_AoryzaeRIB40_AnnotatedCDSs.fasta", "r") as f:
        sequences = [line.strip() for line in f if not line.startswith(">")]
    num_sequences = len(sequences)
    if num_sequences >= 100:
        sequences = sequences[:100]
    up_homology_arm = "AGTCTCGTGTACCGCAGTGGGAGGGCGTGGTGAAGCGAGCTCGGAGGGCGGCGCGCCTCGC"
    down_homology_arm = "CGACTAGTGTCCTGCTGGAGTCGTTCTTGTAGAAAGGAGAAGCCTGCGGCCGCTCGAGCC"
    nucleotide_homology = 30
    generate_primers(sequences, up_homology_arm, down_homology_arm, nucleotide_homology)