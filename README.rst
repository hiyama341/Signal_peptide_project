Signal Peptides
----------------
A signal peptide is a short amino acid sequence(15-30 amino acids) 
that is (typically)present at the N-terminus (beginning) of a protein. 
Its function is to direct the protein to its proper location
within a cell or to be secreted outside of the cell.[1]

Identification
--------------
One widely used prediction algorithm is called signalP, which
uses artificial neural networks to identify signal peptides
based on their amino acid sequence and other features.[2]
Which will also be used in this work. 



Signal peptide architecture
---------------------------
The signal peptide itself can be divided into three regions: the n-region, h-region, and c-region.

.. image:: https://github.com/hiyama341/Signal_peptide_project/blob/cfec5ea8a8000e267c50e2d670d8c413e11e5b33/pictures/Eukaryotic_SP_architecture.png
  :width: 800
  :alt: SP architecture


**n-region**: Determines the orientation of the signal peptide in the protein-conducting channel, Sec61.[1]


**h-region**: The length of this region appears to be a key determinant. Despite the high variability in amino 
acid sequences among signal peptides, the h-region's low sequence conservation suggests it interacts more with 
the SPC's lipid environment than with its transmembrane helices.


**c-region**: Typically 3–7 amino acids long and contains two crucial positions relative to the scissile
bond (−1 and −3), which need to be occupied by small, non-charged residues. 
This is often referred to as the "Ala-X-Ala" rule due to the common presence of alanine (Ala) at these positions. 
Despite the absence of a strong consensus sequence, signal sequences are recognized by 
Signal Peptidase I (SPase I) with high fidelity due to this motif.
However, eukaryotic signal peptides are more flexible, frequently accepting 
Gly, Ser, Thr, and Cys at the −1 position and Ile, Leu, Val, Ser, and Thr at the −3 position almost as often as Ala[3].


Signal Peptidase Complex (SPC) and the interaction with signal peptides
-----------------------------------------------------------------------
Cleavage of signal peptides is performed by the Signal Peptidase Complex (SPC). 
Though our understanding of SPC largely comes from studies in humans, it is likely
that similar mechanisms exist in other eukaryotes[1].

The SPC utilizes a serine protease to cleave the peptide, with a SER-HIS-ASP catalytic triad.

Below you will find a mechanistic understanding of the SPC and their interaction with singal peptides described by Liaci et al. 2021.


.. image:: https://github.com/hiyama341/Signal_peptide_project/blob/81997579cd1d9b3d804a1f55a4fefe1c05291a1a/pictures/signal_peptide_recognition.png
  :width: 800
  :alt: SP architecture

Goal
----
The goal of this project is to engineer signal peptides for 
the secretion of proteins in Aspergillus oryzae using genetic engineering and machine learning approaches, such as SignalP.

Workflow
--------
1. Identify potential signal peptides in Aspergillus oryzae.
2. Perform secretomics analysis and cross-reference with identified signal peptides to filter out and confirm actual signal peptides.
3. Under controlled conditions (i.e., using the same promoter, terminator, and RFP readout), incorporate the top 10 most promising signal peptides.
4. Scramble the last three positions in the best-performing signal peptide, resulting in a theoretical library of about 8,000 mutations. This step is designed to gain a better understanding of how positional changes affect the cleavage activity of the signal peptidase.
5. Apply high-throughput technologies to design, build, and test/screen the strains created from the mutation library.
6. Implement machine learning frameworks to discern the features that constitute the most effective c-region for a signal peptide.

The same methodology will be employed to optimize the other regions of the signal peptide, aiming for an overall enhanced signal peptide.


Colab notebooks
---------------
Below you will find our colab notebooks that describe all the work we made. 


References
----------

1. Liaci et al,. 2021. https://www.sciencedirect.com/science/article/pii/S1097276521006006?via%3Dihub
2. Teufel et al,.2022. https://www.nature.com/articles/s41587-021-01156-3
3. Tuteja,. 2005. https://www.sciencedirect.com/science/article/abs/pii/S000398610500305X