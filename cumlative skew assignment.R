setwd("/Users/swasha/Downloads")
fasta<-read.delim("E coli genome.txt")
fastaDF<-data.frame(fasta)
totalGs<-matrix(nrow=dim(fastaDF)[1], ncol=1)
totalCs<-matrix(nrow=dim(fastaDF)[1], ncol=1)
fastaDF[,1] <- as.character(fastaDF[,1])

for(i in 1:nrow(fastaDF)){  
  genome_split <- strsplit(fastaDF[ i,1], "")[[1]]
  totalGs[i] <- 0
  totalCs[i] <- 0
  for (genome in genome_split) {
    
    if(genome == 'G'){      
      totalGs[i] <- totalGs[i] + 1
    }
    else if(genome == 'C'){
      totalCs[i] <- totalCs[i] + 1 
    }
  }
}

GminusC<-totalGs-totalCs
GplusC<-totalGs+totalCs

GCskew<-GminusC/GplusC
GCskewDF<-data.frame(GCskew)
attach(GCskewDF)
dim(GCskewDF)
w_divided_c<-80/4639675
GCskewDF$WdivC<-seq(from=0.00001724, to=0.00001724*58021, by=0.00001724)
GCskewDF$wc_percent<-GCskewDF$WdivC*100 #X-axis
GCskewDF$GCmultwc<-GCskew*0.00001724    
GCskewDF$cumulative<-0
GCskewDF$cumulative[1]<-GCskewDF$GCmultwc[1]
for(k in 2:nrow(GCskewDF)){
  GCskewDF$cumulative[k]<-GCskewDF$GCmultwc[k]+GCskewDF$cumulative[k-1]
}   #cumulative=Y-axis
write.csv(GCskewDF, file="cumulativeGCskew.csv")
