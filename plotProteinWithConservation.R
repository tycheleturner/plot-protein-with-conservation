#!/bin/R
#Tychele N. Turner
#Laboratory of Aravinda Chakravarti, Ph.D.
#Protein Plotting Script with Conservation
#Programming Language: R
#Updated 06/14/2013

#Description: This script takes mutation information at the protein level and plots out the mutation above the schematic of the protein. It also plots the domains. This version can also add a track for conservation. If you want to use the conservation track the seqinr package will need to be installed in your R instance: install.packages("seqinr"). There are now 2 implementations of this script controlled by the additionalOptions argument. If user specifies yes the other options must be filled out as well. See usage for no (Usage If Without Extra Options) and yes (Usage With Extra Options) answers to the additionalOptions argument.

#NOTE: All files should be referring to the same isoform of the protein. This is imperative for drawing the plot correctly.

#Package requirements: To use the conservation track the seqinr package will need to be installed in your R instance: install.packages("seqinr")

#Required files:
##Mutation file: tab-delimited file containing 5 columns (ProteinId, GeneName, ProteinPositionOfMutation, ReferenceAminoAcid, AlternateAminoAcid) NO HEADER FOR NEEDED FOR THIS FILE
##Protein architecture file: tab-delimited file containing 3 columns (architecture_name, start_site, end_site). This file NEEDS the header and it is the same as what was previously written. This information can be downloaded from the HPRD (http://hprd.org/). Although the most recent files are quite old so looking in the web browser you can get much more up to date information.
##Post-translational modification file: This is a tab-delimited file with only one column and that is the site. This file NEEDS a header and is as previously written.
##Alignment file: This is an aligned multiple sequence alignment fasta file such as that produced by MUSCLE (http://www.ebi.ac.uk/Tools/msa/muscle/). 

#Usage If Without Extra Options:
## R --slave --vanilla < plotProteinWithConservation.R mutationFile proteinArchitectureFile postTranslationalModificationFile alignmentFile referenceSequencePositionInFile nameOfYourQuery tickSize additionalOptions

#Example for usage without extra options:
## R --slave --vanilla < plotProteinWithConservation.R psen1_mutation_file.txt psen1_architecture_file.txt psen1_post_translation_file.txt muscle-I20130227-165316-0600-58424624-pg.fasta 4 Test 20 no

#Usage With Extra Options:
## R --slave --vanilla < plotProteinWithConservation.R mutationFile proteinArchitectureFile postTranslationalModificationFile alignmentFile referenceSequencePositionInFile nameOfYourQuery tickSize additionalOptions showLabels showConservationScore showReferenceSequence showGridlinesAtTicks zoomIn zoomStart zoomEnd wantSecondMutationFileForPlot secondMutationFileForPlot nameOfYourSecondQuery

#Example for usage with extra options:
## R --slave --vanilla < plotProteinWithConservation.R psen1_mutation_file.txt psen1_architecture_file.txt psen1_post_translation_file.txt muscle-I20130227-165316-0600-58424624-pg.fasta 4 Test 2 yes yes yes yes yes yes 110 140 yes psen1_mutation_file.txt Duplicate

#Arguments:
argv <- function(x){
    args <- commandArgs()
    return(args[x])
}

mutationFile <- argv(4) #This is the mutation file
proteinArchitectureFile <- argv(5) #This is the protein architecture file
postTranslationalModificationFile <- argv(6) #This is the post-translation modification file
alignmentFile <- argv(7) #Multiple sequence alignment file
referenceSequencePositionInFile <- argv(8) #This is just the sequence number of your reference in the aligned fasta; for example if human is your reference and its the 5th sequence in the aligned file just say 5
nameOfYourQuery <- argv(9) #Here you can put whatever name you want to show up in the plot
tickSize <- as.numeric(argv(10)) #This is the size of the x-axis tick spacing
additionalOptions <- argv(11) #This requires a yes/no answer, if yes the additional arguments have to be added (see usage), if no, this is the last required argument

if(additionalOptions == "yes"){
	showLabels <- argv(12) #show labels above the mutations
	showConservationScore <- argv(13) #show the actual conservation score 
	showReferenceSequence <- argv(14) #show the sequence of the user-specified reference
	showGridlinesAtTicks <- argv(15) #show full gridlines at each tick mark
	zoomIn <- argv(16) #zoom in to a specific place in the plot (yes/no)
	zoomStart <- argv(17) #zoom start, if zoomIn was no, any value can be added here
	zoomEnd <- argv(18) #zoom end, if zoomIn was no, any value can be added here
	wantSecondMutationFileForPlot <- argv(19) #option to plot a second mutation file
	if(wantSecondMutationFileForPlot == "yes"){ #if you chose yes for plotting a second mutation file the following arguments must be answered
		secondMutationFileForPlot <- argv(20) #file name for the second mutation file
		nameOfYourSecondQuery <- argv(21) #name you want for the second mutation set
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
		if((names(tab[[i]][j])) == a[,5][i])
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
