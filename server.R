#!/bin/R
#Tychele N. Turner
#Laboratory of Aravinda Chakravarti, Ph.D.
#Protein Plotting Script with Conservation Shiny SERVER Script
#Programming Language: R
#Updated 06/14/2013

#Description: This script is the server.R script required for the Shiny application of Plot Protein With Conservation: Visualization of Mutations

library(shiny)
library("seqinr")

# Define server
shinyServer(function(input, output) {

  output$plot <- renderPlot({
    
########mutation file########
	if (is.null(input$mutationFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    var <- read.table(input$mutationFile$datapath, sep="\t")
########mutation file########


########second mutation file########
    if(input$second){
		if (is.null(input$secondmutationFile)) {
    	  # User has not uploaded a file yet
      	return(NULL)
    	}
    
		#cat("input$file=",input$file)
   		 print(paste("input$file$datapath=",input$file$datapath))    
    
    	###dFuserData <- read.csv(input$file$datapath)
	
	    var2 <- read.table(input$secondmutationFile$datapath, sep="\t")
}

########second mutation file########



		
########protein architecture file########
	if (is.null(input$proteinArchitectureFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    pa <- read.table(input$proteinArchitectureFile$datapath, sep="\t", header=TRUE)
########protein architecture file########
	
########post translational modification file########
	if (is.null(input$postTranslationalModificationFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    pt <- read.table(input$postTranslationalModificationFile$datapath, sep="\t", header=TRUE)
########post translational modification file########
	
########conservation file########
	if (is.null(input$alignmentFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    a <- read.fasta(input$alignmentFile$datapath)
########conservation file########	
	
	#CONSERVATION#
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
	referenceSeq <- tdf[which(tdf[,as.numeric(input$referenceSequencePositionInFile)] != "-"),]
	referenceSeq <- as.data.frame(referenceSeq)

	counter <- rep(0, nrow(referenceSeq))
	a <- data.frame(lapply(referenceSeq, as.character), stringsAsFactors=FALSE)
	for(i in 1:nrow(a)){
		a[i,"consensus"] <- paste(as.character(a[i,]), collapse="")
		a[i,"consensus"] <- strtrim(a[i,"consensus"], width=(ncol(a)-1))
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

	
	layout(matrix(c(1,2),nrow=1), widths=c(1.5,4.15))
#	par(oma=c(4, 0, 4, 0), mar=c(5, 0, 4, 0) + 0.1)

	xlimRegion <- c(0, length(scorePlot))

        if(input$zoom) {

                xlimRegion <- c(as.numeric(input$zoomStart), as.numeric(input$zoomEnd))
}

	#stable legend
	plot((-30:-15), rep(-1, 16), col="purple3", type="l", ann=FALSE, bty="n", xaxt="n", yaxt="n", xlim=c(-160, -15), ylim=c(1,-5.5))
		lines((-30:-15), rep(0, 16), col="purple3")
		lines((-30:-15), rep(-0.5, 16), col="purple3")
		text(-100, -0.5, "Conservation", col="purple3", cex=0.9, font=2)
		text(-45,-1, "1", col="purple3", cex=0.9)
		text(-45,-0.5, "0.5", col="purple3", cex=0.9)
		text(-45,0, "0", col="purple3", cex=0.9)
	
	#query text
	text(-100,-2.5,input$nameOfQuery, col="blue", cex=0.9, font=2)

	if(input$second){
		#control text
		#query text
		text(-100,-3.5,input$nameOfQuery2, col="blue", cex=0.9, font=2)
	}

	if(input$showReferenceSequence){
		#reference text
		text(-100,-4.5,"Reference", col="black", cex=0.9, font=2)
	}

	if(input$showConservationScore){
		#score text
		text(-100,0.5,"Score", col="purple3", cex=0.9, font=2)
	}

	plot((1:length(scorePlot)), rep(-2, length(scorePlot)), type="l", lwd=5, main=paste("Amino Acid Changes in", " ", as.character(var[1,2]), " ", "(", as.character(var[1,1]), ")", sep=""), xlab="Amino Acid Position", ylab="", xlim=xlimRegion, ylim=c(1,-5.5), cex.lab=0.9, cex.main=1, yaxt="n", bty="n", font=2, xaxt="n")

	
	#legend
	legend("topleft", legend=c("Protein Domain", "Post-Translational Modification"), fill=c("lightseagreen", "deeppink"), box.col="white", bg="white", cex=1)

	#conservation
	lines(scorePlot, col="purple3")



	#Plot mutations
	points(var[,3], rep(-2.5, length(var[,3])), pch=19, col="blue", cex=0.7)

	#Label mutations
	if (input$labels) {
		for(i in 1:nrow(var)){
			text(var[i,3], rep(-2.6, length(var[i,3])), paste(as.character(var[i,4]), as.character(var[i,3]), as.character(var[i,5]), sep=""), col="blue", cex=0.7, srt=90, adj = 0)
		}
	}


        if(input$second) {
				#Plot mutations
				points(var2[,3], rep(-3.5, length(var2[,3])), pch=19, col="blue", cex=0.7)
				
						#Label mutations
						if (input$labels) {
							for(i in 1:nrow(var2)){
							text(var2[i,3], rep(-3.6, length(var2[i,3])), paste(as.character(var2[i,4]), as.character(var2[i,3]), as.character(var2[i,5]), sep=""), col="blue", cex=0.7, srt=90, adj = 0)
							}
					}
	}


	#labels
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
	
	#Sequence of Interest 
	seqForPlot <- seq[[as.numeric(input$referenceSequencePositionInFile)]][which(seq[[as.numeric(input$referenceSequencePositionInFile)]] != "-")]
	
	if(input$showConservationScore){
		rect(0, 0.3, length(scorePlot), 0.7, col="white", border=NA)
		for(i in 1:length(seqForPlot)){
		text(i,0.5, toupper(abs(round(scorePlot[i],1))), font=2, cex=0.8, srt=90, col="purple3")
		}
	}
			
	if(input$showReferenceSequence){
		rect(0, -4.35, length(scorePlot), -4.65, col="white", border=NA)
		for(i in 1:length(seqForPlot)){
			text(i,-4.5, toupper(seqForPlot[i]), font=2, cex=1)
		}
	
	}
	
	
	
	ticks=seq(0,length(scorePlot), by=as.numeric(input$tickSize))
	#Specify the ticks and grids you want	
	axis(side = 1, at = ticks, las=3)

	if(input$showGridsAtTicks){
		#grid
		for(i in 1:length(seqForPlot)){
			abline(v=ticks[i], lty=3, lwd=0.5)	
		}
	}
		
    })
  #end of plot
  
  
  #called when user uploads a file
   output$table <- renderTable({
########mutation file########
	if (is.null(input$mutationFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    var <- read.table(input$mutationFile$datapath, sep="\t")
########mutation file########
		
########protein architecture file########
	if (is.null(input$proteinArchitectureFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    pa <- read.table(input$proteinArchitectureFile$datapath, sep="\t", header=TRUE)
########protein architecture file########
	
########post translational modification file########
	if (is.null(input$postTranslationalModificationFile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
	#cat("input$file=",input$file)
    print(paste("input$file$datapath=",input$file$datapath))    
    
    ###dFuserData <- read.csv(input$file$datapath)
    pt <- read.table(input$postTranslationalModificationFile$datapath, sep="\t", header=TRUE)
########post translational modification file########
	
	
	#Architecture
	rpa <- pa[order(pa$start_site),]
	for(i in 1:nrow(rpa)){
		rpa$domain_name[i] <- paste(rpa$architecture_name[i], " (", "Domain Number",  i, ") ", sep="")
	}

	var$name <- "Not in Domain"
	for(i in 1:nrow(var)){
		for(j in 1:nrow(rpa)){
			if(as.numeric(var[,3][i]) >= as.numeric(rpa$start_site[j]) & (as.numeric(var[,3][i]) <= as.numeric(rpa$end_site[j])))
			var$name[i] <- rpa$domain_name[j]
		}
	}

	colnames(var) <- c("proteinID", "geneName", "aminoAcidPosition", "refAA", "altAA", "domain")

	#Post Translational Modification
	rpt <- pt[order(pt$site),]

	var$postTranslationSite <- "No"

	for(i in 1:nrow(var)){
		for(j in 1:length(rpt)){
			if(as.numeric(var[,3][i]) == as.numeric(rpt[j]))
			var$postTranslationSite[i] <- "Yes"
		}
	}

	data.frame(var)
	
	
	  
  }) # end of filetable 

 })

