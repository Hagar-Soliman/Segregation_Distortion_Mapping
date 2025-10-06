IM62<-read.delim("/home/hks25/palmer_scratch/libs/aligned/parents/IM62.bam.txt",header=FALSE)
IMPO<-read.delim("/home/hks25/palmer_scratch/libs/aligned/parents/IMPO.bam.txt",header=FALSE)

newgeno <- IM62

count = 0
count2 = 0

# IMPO
for(i in 1:nrow(newgeno)) {
  if(IMPO[i,3] > 0 & IMPO[i,4] > 0) {
    count = count + 1
  }
}

# IM62
for(i in 1:nrow(newgeno)) {
  if(IM62[i,3] > 0 & IM62[i,4] > 0) {
    count2 = count2 + 1
  }
}

count3 = 0
for(i in 1:nrow(newgeno)) {
  if(IM62[i,3] > 0 & IM62[i,4] > 0 & IMPO[i,3] > 0 & IMPO[i,4] > 0) {
    count3 = count3 + 1
  }
}

parents <- matrix(nrow=nrow(IM62), ncol=2)
for(i in 1:nrow(parents)) {
  if(IMPO[i,3] == 0 & IMPO[i,4] != 0) {
    parents[i,1] <- "ALT"
  } else if(IMPO[i,3] != 0 & IMPO[i,4] == 0) {
    parents[i,1] <- "REF"
  } else if(IMPO[i,3] == 0 & IMPO[i,4] == 0) {
    parents[i,1] <- 0
  } else if(IMPO[i,3] != 0 & IMPO[i,4] != 0) {
    parents[i,1] <- "HET"
  }
}

for(i in 1:nrow(parents)) {
  if(IM62[i,3] == 0 & IM62[i,4] != 0) {
    parents[i,2] <- "ALT"
  } else if(IM62[i,3] != 0 & IM62[i,4] == 0) {
    parents[i,2] <- "REF"
  } else if(IM62[i,3] == 0 & IM62[i,4] == 0) {
    parents[i,2] <- 0
  } else if(IM62[i,3] != 0 & IM62[i,4] != 0) {
    parents[i,2] <- "HET"
  }
}

sites <- IMPO[,1:2]
pars <- cbind(sites, parents)
colnames(pars) <- c("Scaffold", "pos", "IMPO", "IM62")
write.table(pars, file="parental_ancestry.tt", row.names=FALSE, quote=FALSE, sep="\t")

count4 = 0
for(i in 1:nrow(pars)) {
  if(pars[i, "IM62"] == pars[i, "IMPO"] & pars[i, "IMPO"] != 0) {
    count4 = count4 + 1
  }
}

count4 = 0
for(i in 1:nrow(pars)) {
  if(pars[i, "IMPO"] == "HET" & pars[i, "IMPO"] == 0) {
    count4 = count4 + 1
  }
}

path <- "/home/hks25/palmer_scratch/libs/aligned/1A_refalt/"
file.names <- dir(path, pattern=".txt")

setwd(path)
output.path <- "/home/hks25/palmer_scratch/libs/aligned/1A_refalt/"

library(parallel) 

# Define a function to process a single file
process_file <- function(file_name) {
  print(paste("Processing file:", file_name))
  indiv <- read.delim(file_name, header=FALSE)
  geno <- data.frame(matrix(nrow=nrow(sites), ncol=2))
  
  for (i in 1:nrow(sites)) {
  if (i %% 1000 == 0) {  # Print progress every 1000 loci
    print(paste("  Processed", i, "loci out of", nrow(sites), "for file:", file_name))
  }
  
  site_chrom <- sites[i, 1]
  site_pos <- sites[i, 2]
  matching_row <- which(indiv[, 1] == site_chrom & indiv[, 2] == site_pos)
  
  if (length(matching_row) == 0) {
    next
  }
    
    indiv_row <- matching_row[1]
    if (indiv[indiv_row, 3] == 0 & indiv[indiv_row, 4] == 0) {
      geno[i, 1:2] <- c(0, 0)
    } else if (pars[i, "IM62"] == pars[i, "IMPO"]) {
      geno[i, 1:2] <- c(0, 0)
    } else if (pars[i, "IMPO"] == "HET" | pars[i, "IM62"] == "HET") {
      geno[i, 1:2] <- c(0, 0)
    } else if (pars[i, "IM62"] == 0 | pars[i, "IMPO"] == 0) {
      geno[i, 1:2] <- c(0, 0)
    } else if (pars[i, "IMPO"] == "REF" & pars[i, "IM62"] == "ALT") {
      geno[i, 1:2] <- indiv[indiv_row, 3:4]
    } else if (pars[i, "IMPO"] == "ALT" & pars[i, "IM62"] == "REF") {
      geno[i, 1] <- indiv[indiv_row, 4]
      geno[i, 2] <- indiv[indiv_row, 3]
    } else {
      print(indiv[indiv_row, 3:4])
      print(pars[, 3:4])
    }
  }
  
  data <- cbind(sites, geno)
  output_file <- paste0(output.path, basename(file_name))
  write.table(data, file=output_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

# Use mclapply to process files in parallel
num_cores <- 8  
print(paste("Using", num_cores, "cores for parallel processing"))

# Use mclapply to process files in parallel
mclapply(file.names, process_file, mc.cores=num_cores)
