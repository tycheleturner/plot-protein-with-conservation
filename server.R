library(shiny)
library("seqinr")

# Define server logic for random distribution application
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
	#x is the input data, y is rpt, z is rpa from HPRD
	#pdf(paste(as.character(var[1,2]), "_protein_plot.pdf", sep=""), height=7.5, width=10)
	#par(oma=c(8, 1.2, 8, 1.2))

	xlimRegion <- c(-145, length(scorePlot))

        if(input$zoom) {

                xlimRegion <- c(as.numeric(input$zoomStart), as.numeric(input$zoomEnd))
}

	plot((1:length(scorePlot)), rep(-2, length(scorePlot)), type="l", lwd=5, main=paste("Amino Acid Changes in", " ", as.character(var[1,2]), " ", "(", as.character(var[1,1]), ")", sep=""), xlab="Amino Acid Position", ylab="", xlim=xlimRegion, ylim=c(0,-4), cex.lab=0.9, cex.main=1, yaxt="n")

	#Plot mutations
	points(var[,3], rep(-2.5, length(var[,3])), pch=19, col="blue", cex=0.7)

	#Label mutations
	if (input$labels) {
		for(i in 1:nrow(var)){
			text(var[i,3], rep(-2.7, length(var[i,3])), paste(as.character(var[i,4]), as.character(var[i,3]), as.character(var[i,5]), sep=""), col="blue", cex=0.7, srt=90, adj = 0)
		}
	}

	#labels
	text(-100,-2.5,input$nameOfQuery, col="blue", cex=0.9)
	for(i in 1:length(pt$site)){
		segments(as.numeric(pt$site[i]), -2, as.numeric(pt$site[i]), -2.25, lwd=2, col="black")
		points(as.numeric(pt$site[i]), -2.25, pch=19, col="deeppink", cex=0.7)
	}

	for(i in 1:length(pa$start_site)){
#		lines(c(as.numeric(pa$start_site[i]):as.numeric(pa$end_site[i])), rep(-2,(as.numeric(pa$end_site[i])- as.numeric(pa$start_site[i]))+1),lwd=10, col=c("purple"))
	rect(as.numeric(pa$start_site[i]), -2.1, as.numeric(pa$end_site[i]), -1.9, col="lightseagreen")
	}
	for(i in 1:length(pa$architecture_name)){
		text(median(c(as.numeric(pa$start_site[i]), as.numeric(pa$end_site[i]))), -1.80, pa$architecture_name[i], cex=0.6)
	}
	lines(scorePlot, col="purple3")
	lines((-150:-40), rep(-1, 111), col="purple3")
	lines((-150:-40), rep(0, 111), col="purple3")
	text(-95, -0.5, "Conservation", col="purple3", cex=0.9)
	text(-20,-1, "1", col="purple3", cex=0.9)
	text(-20,0, "0", col="purple3", cex=0.9)
	legend("topright", c("Protein Domain", "Post-Translational Modification"), fill=c("lightseagreen", "deeppink"), cex=0.7)
	#dev.off()


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

