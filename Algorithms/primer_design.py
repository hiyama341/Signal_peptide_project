"""Function used to generate forward and reverse primers from positive strand nucleotide sequences"""
import pandas as pd
from typing import List
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def generate_primers(nucleotide_sequences: List[SeqRecord or str], up_homology_arm: str = None,
                                  down_homology_arm: str = None, nucleotide_homology: int = None) -> pd.DataFrame:
    """
    This function generates sequencing-ready forward and reverse primers by processing a list of String or SeqRecords nucleotide sequences 
    defined in 5´ to 3´ direction from the positive (top) strand of the nucleotide sequence to be expressed. 
        The process consists of several steps including:
            1. Addition of up and down homology/repair arms to the nucleotide sequence to create homology with the bakcbone sequenece
            2. Division of the nucleotide sequence into equal halves while preserving codon structure
            3. Sharing defined number of sequences between the nucleotide sequence halves to create matching homology between them
            4. Construction of reverse complement ( 5´ to 3´ direction) sequence of the forward primer on the negative (bottom) strand 
            5. Calculation of the forward and reverse primer lengths
            5. Construction of a dataframe to summarize forward and reverse primers and their lengths

    Parameters:
    -----------
        nucleotide_sequenceuences: str or list of Bio.SeqReccord.SeqRecord
            A list of String or SeqRecords nucleotide sequences defined in 5´ to 3´ direction from the positive (top) strand
        
        up_homology_arm: str
            A string representing nucleotide sequence to be added to the beginning of the forward primer to create homology with the bakcbone sequenece
        
        down_homology_arm: str
            A string representing nucleotide sequence to be added to the end of the reverse primer to create homology with the bakcbone sequence
        
        nucleotide_homology: int
            An integer representing the number of nucleotides that should be shared between the halves of the nucleotide sequence to create homology between them
        

    Returns:
    --------
        A pandas dataframe with forward and reverse primer sequnces and their respective primer lengths.
    """
    # Extract the sequencecs if they are provided as SeqRecord 
    extracted_sequences = []
    
    for sequence in nucleotide_sequences:
        if isinstance(sequence, SeqRecord):
            nucleotide_sequence = str(sequence.seq)
        elif isinstance(sequence, str):
            nucleotide_sequence = sequence
        else:
            raise TypeError(f"Expected SeqRecord, Seq object or string, but got {type(sequence)}")

        # Check if the sequence contains only nucleotides
        if set(nucleotide_sequence) - set('agtACGT'):
            raise ValueError(f"Sequence {sequence} at index {nucleotide_sequences.index(sequence)} contains non-nucleotide characters.")
        
        # Add up and down homology arms to the nuclceotide sequences (if stated)
        if up_homology_arm and down_homology_arm:
            nucleotide_sequence = up_homology_arm + nucleotide_sequence + down_homology_arm
        else:
            nucleotide_sequence = nucleotide_sequence
        
        # Determine the half of the sequence and divide into two while preseving codon structure 
        sequence_length = len(nucleotide_sequence)
        sequence_midpoint = sequence_length // 2  # Find the midpoint of the sequence

        if sequence_midpoint % 3 != 0:  # Check if sequence_midpoint is not divisible by 3
            codon_position = sequence_midpoint % 3  # Find the position of the last complete codon before the midpoint
            if codon_position <= 1:  # Adjust midpoint to conserve codons
                sequence_midpoint -= codon_position
            else:
                sequence_midpoint += (3 - codon_position)

        up_sequence = nucleotide_sequence[:sequence_midpoint]
        down_sequence = nucleotide_sequence[sequence_midpoint:]
        
        # Create homology arms between halves (up and down) nucleotide sequences and emphasize the homology match by lowercase letters    
        up_seq_homology = up_sequence[-nucleotide_homology:].lower()
        down_seq_homology = down_sequence[:nucleotide_homology].lower()

        up_primer = up_sequence + down_seq_homology
        down_primer = up_seq_homology + down_sequence

        # Define the complementary translation table for nucleotides
        nucleotide_complements = str.maketrans("ATCGatcg", "TAGCtagc")
        complementary_up_primer = up_primer.translate(nucleotide_complements) # Generate complementary strand of the nucleotide sequence (3´ to 5´ on negative strand)
        reverse_complementary_up_primer = complementary_up_primer[::-1] # Take the reverse complement of the complementary strand (write in the direction of 5´ to 3´ on negative strand)
       
        # Emphasize the entire homology section between nucleotide sequences by making them lowercase 
        forward_primer = reverse_complementary_up_primer[:2*nucleotide_homology].lower() + reverse_complementary_up_primer[2*nucleotide_homology:]
        reverse_primer = down_primer[:2*nucleotide_homology].lower() + down_primer[2*nucleotide_homology:]

        extracted_sequences.append((forward_primer, len(forward_primer), reverse_primer, len(reverse_primer)))

    df = pd.DataFrame(extracted_sequences, columns=['Forward primer', 'Forward primer length', 'Reverse primer', 'Reverse primer length'])
    return df