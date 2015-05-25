#!/bin/R
#Tychele N. Turner
#Laboratory of Aravinda Chakravarti, Ph.D.
#Protein Plotting Script with Conservation
#Programming Language: R
#Updated 05/25/2015

#Description: This script takes mutation information at the protein level and plots out the mutation above the schematic of the protein. It also plots the domains. This version can also add a track for conservation. If you want to use the conservation track the seqinr package will need to be installed in your R instance: install.packages("seqinr"). There are now 2 implementations of this script controlled by the additionalOptions argument. If user specifies yes the other options must be filled out as well. See usage for no (Usage If Without Extra Options) and yes (Usage With Extra Options) answers to the additionalOptions argument.

#NOTE: All files should be referring to the same isoform of the protein. This is imperative for drawing the plot correctly.

#Package requirements: To use the conservation track the seqinr package will need to be installed in your R instance: install.packages("seqinr")

#Usage
##Basic
###Rscript plotProteinWithConservation.R -m psen1_mutation_file.txt -a psen1_architecture_file.txt -p psen1_post_translation_file.txt -f muscle-I20130227-165316-0600-58424624-pg.fasta -r 4

##Advanced
###Rscript plotProteinWithConservation.R -m psen1_mutation_file.txt -a psen1_architecture_file.txt -p psen1_post_translation_file.txt -f muscle-I20130227-165316-0600-58424624-pg.fasta -r 4 -n Disease -t 25 -v yes -s yes -d yes -e yes -j yes -z yes -b 50 -c 100 -q yes -u psen1_mutation_file.txt -y Disease2

#will have to run the install for optparse if never installed before
#install.packages("optparse")

library("optparse")

option_list <- list(
    make_option(c('-m', '--mutations'), action='store', type='character', default='mutationFile.txt', help='This is the mutation file. It should be a tab-delimited file containing 5 columns (ProteinId, GeneName, ProteinPositionOfMutation, ReferenceAminoAcid, AlternateAminoAcid) NO HEADER FOR NEEDED FOR THIS FILE. (REQUIRED)'),
    make_option(c('-a', '--architecture'), action='store', type='character', default='architectureFile.txt', help='This is the protein architecture file. It should be a tab-delimited file containing 3 columns (architecture_name, start_site, end_site). This file NEEDS the header and it is the same as what was previously written. This information can be downloaded from the HPRD (http://hprd.org/). Although the most recent files are quite old so looking in the web browser you can get much more up to date information. (REQUIRED)'),
    make_option(c('-p', '--posttranslational'), action='store', type='character', default='posttranslationalFile.txt', help='This is the protein post-translational modification file. This is a tab-delimited file with only one column and that is the site. This file NEEDS a header and is as previously written (site). (REQUIRED)'),
    make_option(c('-f', '--fastaalignmentfile'), action='store', type='character', default='test.fasta', help='This is an aligned multiple sequence alignment fasta file such as that produced by MUSCLE (http://www.ebi.ac.uk/Tools/msa/muscle/). (REQUIRED)'),
    make_option(c('-r', '--referencesequenceposition'), action='store', type='numeric', default=1, help='Reference sequence position in the fasta alignment file (REQUIRED)'),
    make_option(c('-n', '--name'), action='store', type='character', default='Test', help='Name of your query. Default is Test'),
    make_option(c('-t', '--ticksize'), action='store', type='numeric', default=10, help='Size of ticks on x-axis. Default is 10'),
    make_option(c('-v', '--advancedoptions'), action='store', type='character', default='no', help='Whether or not to look at advanced options including showing labels above mutations, showing actual conservation score, showing sequence of user-specified reference, showing full gridlines at each tick mark, zooming in to different parts of the plot, zoom starts, zoom ends, and showing a second mutation set. Default is no'),
    make_option(c('-s', '--showlabels'), action='store', type='character', default='no', help='Option to show labels. Default is no'),
    make_option(c('-d', '--showconservationscore'), action='store', type='character', default='no', help='Option to show conservation score. Default is no'),
    make_option(c('-e', '--showreferencesequence'), action='store', type='character', default='no', help='Option to show reference sequence of user-defined reference. Default is no'),
    make_option(c('-j', '--showgridlines'), action='store', type='character', default='no', help='Option to show gridlines at each tick. Default is no'),
    make_option(c('-z', '--zoom'), action='store', type='character', default='no', help='Option to zoom in. Default is no'),
    make_option(c('-b', '--zoomstart'), action='store', type='numeric', default=1, help='Starting number for zoom in. Use if zoom option is set to yes. Default is 1'),
    make_option(c('-c', '--zoomend'), action='store', type='numeric', default=10, help='Ending number for zoom in. Use if zoom option is set to yes. Default is 10'),
    make_option(c('-q', '--secondmutationfile'), action='store', type='character', default='no', help='Option to enter a second mutation file. Default is no'),
    make_option(c('-u', '--secondmutationfilename'), action='store', type='character', default='test2.txt', help='Name of second mutation plot file'),
    make_option(c('-y', '--secondmutationsetname'), action='store', type='character', default='Test2', help='Name of second mutation dataset. Default is Test2')
)
opt <- parse_args(OptionParser(option_list = option_list))

mutationFile <- opt$mutations
proteinArchitectureFile <- opt$architecture
postTranslationalModificationFile <- opt$posttranslational
alignmentFile <- opt$fastaalignmentfile
referenceSequencePositionInFile <- opt$referencesequenceposition
nameOfYourQuery <- opt$name
tickSize <- opt$ticksize
additionalOptions <- opt$advancedoptions

if(additionalOptions == "yes"){
	showLabels <- opt$showlabels
	showConservationScore <- opt$showconservationscore
	showReferenceSequence <- opt$showreferencesequence
	showGridlinesAtTicks <- opt$showgridlines
	zoomIn <- opt$zoom
	zoomStart <- opt$zoomstart
	zoomEnd <- opt$zoomend
	wantSecondMutationFileForPlot <- opt$secondmutationfile
	if(wantSecondMutationFileForPlot == "yes"){
		secondMutationFileForPlot <- opt$secondmutationfilename
		nameOfYourSecondQuery <- opt$secondmutationsetname
		}
	}

####################ANALYSIS####################
library("seqinr")
#Read in the files
var <- read.table(mutationFile, sep="\t")
pa <- read.table(proteinArchitectureFile, sep="\t", header=TRUE)
pt <- read.table(postTranslationalModificationFile, sep="\t", header=TRUE)

if(additionalOptions == "yes"){
	if(wantSecondMutationFileForPlot == "yes"){
		var2 <- read.table(secondMutationFileForPlot, sep="\t")
	}
}

#conservation analysis
a <- read.fasta(alignmentFile)
seq <- list()
for(i in 1:length(a)){
  seq[[i]] <- a[[i]][1:length(a[[i]])]
}
numberOfSeq <- length(seq)

mat <- matrix(0, nrow=length(a), ncol=length(a[[1]]))
for(i in 1:length(seq)){
	mat[i,] <- seq[[i]]
}

df <- as.data.frame(mat)
tdf <- t(df)
referenceSeq <- tdf[which(tdf[,as.numeric(referenceSequencePositionInFile)] != "-"),]
referenceSeq <- as.data.frame(referenceSeq)
write.table(referenceSeq, file="alignment_table", sep="\t", quote=F, row.names=F, col.names=F)

counter <- rep(0, nrow(referenceSeq))
a <- read.table("alignment_table", sep="\t")
a <- data.frame(lapply(a, as.character), stringsAsFactors=FALSE)
for(i in 1:nrow(a)){
	a[i,"consensus"] <- paste(as.character(a[i,]), collapse="")
	}

countBases <- function(string){
 table(strsplit(string, "")[[1]])
}
c <- as.character(a[,"consensus"])
tab <- list()
for(i in 1:length(c)){
	tab[[i]] <- countBases(c[i])
}

score <- rep(0, nrow(a))
for(i in 1:length(tab)){
	for(j in 1:length(tab[[i]])){
#		if((names(tab[[i]][j])) == a[,5][i])
		if((names(tab[[i]][j])) == a[i,][as.numeric(referenceSequencePositionInFile)])
		score[i] <- tab[[i]][j]
	}
}
scorePlot <- -(((score/numberOfSeq)))
system("rm alignment_table")
############PLOTTING#############
#x is the input data, y is rpt, z is rpa from HPRD
pdf(paste(as.character(var[1,2]), "_protein_plot.pdf", sep=""), height=8.5, width=11)
layout(matrix(c(1,2),nrow=1), widths=c(1,3))
par(oma=c(4, 0, 4, 0), mar=c(5, 0, 4, 0) + 0.4)

#stable legend
plot((-30:-15), rep(-1, 16), col="purple3", type="l", ann=FALSE, bty="n", xaxt="n", yaxt="n", xlim=c(-160, -15), ylim=c(1,-5.5))
	lines((-30:-15), rep(0, 16), col="purple3")
	lines((-30:-15), rep(-0.5, 16), col="purple3")
	text(-100, -0.5, "Conservation", col="purple3", cex=0.9, font=2)
	text(-45,-1, "1", col="purple3", cex=0.9)
	text(-45,-0.5, "0.5", col="purple3", cex=0.9)
	text(-45,0, "0", col="purple3", cex=0.9)
	

#query text
text(-100,-2.5,nameOfYourQuery, col="blue", cex=0.9, font=2)

if(additionalOptions == "yes"){
	if(wantSecondMutationFileForPlot == "yes"){
		text(-100,-3.5,nameOfYourSecondQuery, col="blue", cex=0.9, font=2)
		}
}

if(additionalOptions == "yes"){
	if(showReferenceSequence == "yes"){
		#reference text
		text(-100,-4.5,"Reference", col="black", cex=0.9, font=2)
	}
}

if(additionalOptions == "yes"){
	if(showConservationScore == "yes"){
		#score text
		text(-100,0.5,"Score", col="purple3", cex=0.9, font=2)
	}
}


xlimRegion <- c(0, length(scorePlot))

	if(additionalOptions == "yes"){
		if(zoomIn == "yes") {
    	      xlimRegion <- c(as.numeric(zoomStart), as.numeric(zoomEnd))
		}
	}

#main plot
plot((1:length(scorePlot)), rep(-2, length(scorePlot)), type="l", lwd=5, main=paste("Amino Acid Changes in", " ", as.character(var[1,2]), " ", "(", as.character(var[1,1]), ")", sep=""), xlab="Amino Acid Position", ylab="", xlim=xlimRegion, ylim=c(1,-5.5), cex.lab=0.9, cex.main=1, yaxt="n", bty="n", font=2, xaxt="n")

#legend
legend("topleft", legend=c("Protein Domain", "Post-Translational Modification"), fill=c("lightseagreen", "deeppink"), box.col="white", bg="white", cex=1)

#conservation
lines(scorePlot, col="purple3")

#Sequence of Interest 
seqForPlot <- seq[[as.numeric(referenceSequencePositionInFile)]][which(seq[[as.numeric(referenceSequencePositionInFile)]] != "-")]

ticks=seq(0,length(scorePlot), by=tickSize) 
axis(side = 1, at = ticks, las=3)


if(additionalOptions == "yes"){	
	if(showGridlinesAtTicks == "yes"){
		#grid
		for(i in 1:length(seqForPlot)){
		abline(v=ticks[i], lty=3, lwd=0.5, col="lightgray")
		}
	}
}

#Plot mutations
points(var[,3], rep(-2.5, length(var[,3])), pch=19, col="blue", cex=0.7)

if(additionalOptions == "yes"){
	if(showLabels == "yes"){
		#Label mutations
		for(i in 1:nrow(var)){
			text(var[i,3], rep(-2.6, length(var[i,3])), paste(as.character(var[i,4]), as.character(var[i,3]), as.character(var[i,5]), sep=""), col="blue", cex=0.8, srt=90, adj = 0)
			}
		}	
}

if(additionalOptions == "yes"){
	if(wantSecondMutationFileForPlot == "yes"){
		#Plot mutations
		points(var2[,3], rep(-3.5, length(var2[,3])), pch=19, col="blue", cex=0.7)

		if(showLabels == "yes"){
			#Label mutations
			for(i in 1:nrow(var2)){
			text(var2[i,3], rep(-3.6, length(var2[i,3])), paste(as.character(var2[i,4]), as.character(var2[i,3]), as.character(var2[i,5]), sep=""), col="blue", cex=0.8, srt=90, adj = 0)
			}
		}
	}
}

#labels post-translation and protein architecture
for(i in 1:length(pt$site)){
	segments(as.numeric(pt$site[i]), -2, as.numeric(pt$site[i]), -2.32, lwd=2, col="black")
	points(as.numeric(pt$site[i]), -2.32, pch=19, col="deeppink", cex=0.7)
}

for(i in 1:length(pa$start_site)){
	rect(as.numeric(pa$start_site[i]), -2.1, as.numeric(pa$end_site[i]), -1.9, col="lightseagreen")
}
for(i in 1:length(pa$architecture_name)){
	text(median(c(as.numeric(pa$start_site[i]), as.numeric(pa$end_site[i]))), -1.70, pa$architecture_name[i], cex=1)
}

if(additionalOptions == "yes"){
	if(showReferenceSequence == "yes"){
		#specify if you want a reference track
		rect(0, -4.35, length(scorePlot), -4.65, col="white", border=NA)
		for(i in 1:length(seqForPlot)){
			text(i,-4.5, toupper(seqForPlot[i]), font=2, cex=1)
			}
		}
	}
	
if(additionalOptions == "yes"){
	if(showConservationScore == "yes"){

		#specify if you want score
		rect(0, 0.3, length(scorePlot), 0.7, col="white", border=NA)
		for(i in 1:length(seqForPlot)){
			text(i,0.5, toupper(abs(round(scorePlot[i],1))), font=2, cex=0.8, srt=90, col="purple3")
			}
		}
	}

dev.off()
