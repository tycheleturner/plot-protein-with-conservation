plot-protein-with-conservation
==============================

Plot Protein: Visualization of Mutations with Conservation

Author: Tychele N. Turner, Laboratory of Aravinda Chakravarti, Ph.D.

Licenses: GNU General Public License version 3.0 (GPLv3), MIT License

Short Description: Protein Plotting Script with Conservation

Programming Language: R

Readme Date: 04/19/2013

Description: This script takes mutation information at the protein level and plots out the mutation above the schematic of the protein. It also plots the domains. This version can also add a track for conservation. If you want to use the conservation track the seqinr package will need to be installed in your R instance: install.packages("seqinr")

NOTE: All files should be referring to the same isoform of the protein. This is imperative for drawing the plot correctly.

Package requirements: To use the conservation track the seqinr package will need to be installed in your R instance: install.packages("seqinr")

Required files:

*Mutation file: tab-delimited file containing 5 columns (ProteinId, GeneName, ProteinPositionOfMutation, ReferenceAminoAcid, AlternateAminoAcid) NO HEADER FOR NEEDED FOR THIS FILE

*Protein architecture file: tab-delimited file containing 3 columns (architecture_name, start_site, end_site). This file NEEDS the header and it is the same as what was previously written. This information can be downloaded from the HPRD (http://hprd.org/). Although the most recent files are quite old so looking in the web browser you can get much more up to date information.

*Post-translational modification file: This is a tab-delimited file with only one column and that is the site. This file NEEDS a header and is as previously written.

*Alignment file: This is an aligned multiple sequence alignment fasta file such as that produced by MUSCLE (http://www.ebi.ac.uk/Tools/msa/muscle/). 


Usage:

R --slave --vanilla < plotProteinWithConservation.R mutationFile proteinArchitectureFile postTranslationalModificationFile alignmentFile referenceSequencePositionInFile nameOfYourQuery

Example:

R --slave --vanilla < scripts/plotProteinWithConservation.R testFiles/psen1_mutation_file.txt testFiles/psen1_architecture_file.txt testFiles/psen1_post_translation_file.txt testFiles/muscle-I20130227-165316-0600-58424624-pg.fasta 4 Test
