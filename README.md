plot-protein-with-conservation
==============================

Plot Protein: Visualization of Mutations with Conservation

Author: Tychele N. Turner, Laboratory of Aravinda Chakravarti, Ph.D.

Licenses: GNU General Public License version 3.0 (GPLv3), MIT License

Short Description: Protein Plotting Script with Conservation

Programming Language: R

Version: 2.0.0

Readme Date: 05/25/2015

Description: This script takes mutation information at the protein level and plots out the mutation above the schematic of the protein. It also plots the domains. This version can also add a track for conservation. If you want to use the conservation track the seqinr package will need to be installed in your R instance: install.packages("seqinr") There are now 2 implementations of this script controlled by the additionalOptions argument. If user specifies yes the other options must be filled out as well. See usage for no (Usage If Without Extra Options) and yes (Usage With Extra Options) answers to the additionalOptions argument.

NOTE: All files should be referring to the same isoform of the protein. This is imperative for drawing the plot correctly.

Package requirements: To use the conservation track the seqinr package will need to be installed in your R instance: install.packages("seqinr")

Required files:

*Mutation file: tab-delimited file containing 5 columns (ProteinId, GeneName, ProteinPositionOfMutation, ReferenceAminoAcid, AlternateAminoAcid) NO HEADER FOR NEEDED FOR THIS FILE

*Protein architecture file: tab-delimited file containing 3 columns (architecture_name, start_site, end_site). This file NEEDS the header and it is the same as what was previously written. This information can be downloaded from the HPRD (http://hprd.org/). Although the most recent files are quite old so looking in the web browser you can get much more up to date information.

*Post-translational modification file: This is a tab-delimited file with only one column and that is the site. This file NEEDS a header and is as previously written.

*Alignment file: This is an aligned multiple sequence alignment fasta file such as that produced by MUSCLE (http://www.ebi.ac.uk/Tools/msa/muscle/). 


Basic Usage:
==================================================
Rscript plotProteinWithConservation.R -m psen1_mutation_file.txt -a psen1_architecture_file.txt -p psen1_post_translation_file.txt -f muscle-I20130227-165316-0600-58424624-pg.fasta -r 4

Advanced usage:
==================================================
Rscript plotProteinWithConservation.R -m psen1_mutation_file.txt -a psen1_architecture_file.txt -p psen1_post_translation_file.txt -f muscle-I20130227-165316-0600-58424624-pg.fasta -r 4 -n Disease -t 25 -v yes -s yes -d yes -e yes -j yes -z yes -b 50 -c 100 -q yes -u psen1_mutation_file.txt -y Disease2
