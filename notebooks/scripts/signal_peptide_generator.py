#!/usr/bin/env python
# MIT License
# Copyright (c) 2023, Technical University of Denmark (DTU)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from Bio.Seq import Seq
from typing import List
import numpy as np
import h2o.estimators

'''This is a module to work with signal peptides and generate artificial ones'''

def get_signal_peptides_cross_ref_with_genome(list_of_peptides: List[dict], all_proteins: List[SeqRecord]) -> List[SeqRecord]:
    """
    Extracts the protein sequence that corresponds to each predicted signal peptide sequence from the input list
    of peptides and matches the signal peptide to its corresponding protein sequence in the input list of 
    protein sequences.
    
    Parameters
    ----------
    list_of_peptides : list
        A list of dictionaries containing information about predicted signal peptide sequences, including gene name, 
        start and end positions, and the signal peptide likelihood score.
    all_proteins : list
        A list of SeqRecord objects containing protein sequences.
    
    Returns
    -------
    list
        A list of SeqRecord objects that correspond to the input predicted signal peptide sequences, including 
        protein sequence, ID, name, and a description indicating that the sequence corresponds to a predicted 
        signal peptide.
    """
    signal_peptide_seqs = []

    for signal_peptide in list_of_peptides:
        for seqrecord in all_proteins:
            if signal_peptide['gene'] in seqrecord.id:             
                seq = SeqRecord(
                    Seq(seqrecord.seq[signal_peptide['start_pos']:signal_peptide['end_pos']]), 
                    id=seqrecord.id,
                    name=seqrecord.name,
                    description="signal_peptide predicted by signalP")

                signal_peptide_seqs.append(seq)

    return signal_peptide_seqs


def generate_artificial_peptide(list_of_probabilities: np.ndarray, amino_acids: np.ndarray, max_length=22) -> str:
    """
    Generate an artificial peptide based on a list of probabilities and amino acids.
    
    Parameters:
    ----------
    list_of_probabilities : numpy.ndarray
        2-D array of probability of amino acids in the peptide
    amino_acids : numpy.ndarray
        1-D array of amino acids.
        
    Returns:
    -------
    str
        Generated artificial peptide
        
    Notes:
    ------
    The length of the probability array should be same as the length of the peptide.
    """
    out_str = ''
    for i in range(len(list_of_probabilities)):
        # make synthetic signal peptide
        artificial_amino_acid = list(np.random.choice(amino_acids, 1, p=list_of_probabilities[i]))

        if artificial_amino_acid == ['-']: 
            break

        out_str += artificial_amino_acid[0]
    return out_str



def add_dunder_tail(peptide: str, max_lenght: int = 22) -> str:
    '''
    Adds a tail if the length of the given peptide is less than the specified maximum length.

    Parameters:
    -----------
    peptide: str
        The peptide to which a tail needs to be added.
    max_lenght: int, optional (default=22)
        The maximum length of the peptide including the tail.

    Returns:
    --------
    str
        The peptide with added tail if it was shorter than the specified maximum length.
    '''
    # Check if the length of the given peptide is less than the specified maximum length.
    if len(peptide) < max_lenght:
        # Calculate the difference between the length of the given peptide and the specified maximum length.
        difference = max_lenght - len(peptide)
        # Add '-' to the peptide to create the tail.
        sequence = peptide + ('-' * difference)
    else:
        sequence = peptide

    return sequence

def generate_artificial_peptides(list_of_probabilities: np.ndarray, amino_acids: np.ndarray, n_peptides: int, max_len=50) -> pd.DataFrame:
    """
    Generate a dataframe of artificial peptides based on a list of probabilities and amino acids.
    
    Parameters:
    ----------
    list_of_probabilities : numpy.ndarray
        2-D array of probability of amino acids in the peptide
    amino_acids : numpy.ndarray
        1-D array of amino acids.
    n_peptides : int
        Number of peptides to generate
        
    Returns:
    -------
    pd.DataFrame
        Dataframe of generated artificial peptides with 'sequence' as column
        
    Notes:
    ------
    The length of the probability array should be same as the length of the peptide.
    """
    artificial_peptides = []
    lengths = [] 
    for i in range(n_peptides): 
        peptide = generate_artificial_peptide(list_of_probabilities,amino_acids, max_length=max_len)
        if len(peptide) <= max_len:
            peptide_w_tail = add_dunder_tail(peptide, max_lenght = max_len)
        else: 
            continue
        
        # save
        lengths.append(len(peptide))                                     
        artificial_peptides.append(peptide_w_tail)

    df = pd.DataFrame(artificial_peptides, columns =['sequence'])
    df['length'] = lengths
    return df



def split_peptides_sequences(df_signalPP: pd.DataFrame) -> pd.DataFrame:
    '''
    Splits the amino acid sequences into individual amino acids for each position.

    Parameters:
    -----------
    df_signalPP: pandas.DataFrame
        A DataFrame containing the amino acid sequences.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the split amino acid sequences.
    '''
    # Initialize an empty list to store the split sequences.
    peptides_split = []

    # Split each amino acid sequence into individual amino acids.
    for k, v in df_signalPP.iterrows():
        sequence = []
        for seq in v['sequence']:
            sequence.append(seq)
        peptides_split.append(sequence)

    # Convert the list of split sequences into a DataFrame and fill NaN values with '-'.
    new_peptides = pd.DataFrame(peptides_split)
    new_peptides = new_peptides.fillna('-')

    return new_peptides




def signal_peptide_predictor(list_of_probabilities: list, amino_acids: str, n_peptides: int, number_of_iterations: int, best_model: h2o.estimators, 
                             training_column_name:str = 'MM_N_peptide_abundance', max_len_of_signal_peptides = 30) -> pd.DataFrame:
    '''
    Predicts the best signal peptides from a given number of iterations.

    Parameters:
    -----------
    list_of_probabilities: list
        A list of probabilities.
    amino_acids: str
        A string containing the amino acids used to generate the peptides.
    n_peptides: int
        The number of peptides to generate.
    number_of_iterations: int
        The number of iterations to run the predictor.
    best_model: h2o.estimators
        The trained model to use for prediction.
    training_column_name : str
        The name of the column that the model has been trained on.
    max_len_of_signal_peptides : int
        The maximum lenght of the signal peptides you want to be generated 

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the predicted signal peptides.
    '''
    # Initialize an empty DataFrame to store the predicted peptides.
    data = pd.DataFrame()
    
    # Generate and predict peptides for each iteration.
    for i in range(0, number_of_iterations):
        # Generate new peptides.
        new_TO_NATURE_peptides = generate_artificial_peptides(list_of_probabilities, amino_acids, n_peptides=n_peptides, max_len=max_len_of_signal_peptides)
        
        # Split the peptides into sequences.
        new_TO_NATURE_peptides = split_peptides_sequences(new_TO_NATURE_peptides)

        # Convert the DataFrame to an H2OFrame and make the columns categorical.
        df_test = h2o.H2OFrame(pd.concat([new_TO_NATURE_peptides], axis='columns'))
        for column in df_test.columns:
            if column != training_column_name:
                df_test[column] = df_test[column].asfactor()

        # Make predictions on the test data.
        predicted = best_model.predict(df_test).as_data_frame()
        new_TO_NATURE_peptides['predictions'] = predicted['predict'].to_list()

        # Concatenate the new predictions with the existing DataFrame.
        if len(data) == 0:
            data = new_TO_NATURE_peptides.copy()
        else:
            data = pd.concat([data, new_TO_NATURE_peptides], axis=0)
            data = data.sort_values('predictions', ascending=False)
            data = data[0:100]

        # Print the current iteration number every 1000 iterations.
        if i % 1000 == 0:
            print(f"Iteration {i}")
    
    return data