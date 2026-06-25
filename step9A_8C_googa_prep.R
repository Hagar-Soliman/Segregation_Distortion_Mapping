setwd("/home/hks25/palmer_scratch/googaV3/refAlt/8C_refAlt/50kb_windows/")
path="/home/hks25/palmer_scratch/googaV3/refAlt/8C_refAlt/50kb_windows/"

files<-dir(path, pattern="Genostats\\.")
data<-data.frame(matrix(ncol=5,nrow=length(files)))
for(k in 1:length(files)){
  data[k,1:5]<-read.delim(files[k],header=TRUE)
}
names(data)<-c("indiv","AA","AB","BB","NN")


genotype.files<-dir(path, pattern="Genotypes\\.")
countAA<-0
countAB<-0
countBB<-0
countNN<-0

test<-read.delim(dir(path, pattern="Genotypes\\.")[1], header=TRUE)
n_windows<-nrow(test)

site.counts<-data.frame(matrix(nrow=n_windows,ncol=6))
window.names<-test[,2:3]

for(d in 1:n_windows){
  countAA<-0
  countAB<-0
  countBB<-0
  countNN<-0
  for(m in genotype.files){
    temp<-read.delim(m,header=TRUE)
    if(temp[d,4]=="AA"){
      countAA=countAA+1
    } else if(temp[d,4]=="AB"){
      countAB=countAB+1
    } else if(temp[d,4]=="BB"){
      countBB=countBB+1
    } else if(temp[d,4]=="NN"){
      countNN=countNN+1
    }
  }
  site.counts[d,3]<-countAA
  site.counts[d,4]<-countAB
  site.counts[d,5]<-countBB
  site.counts[d,6]<-countNN
}
site.counts[,1:2]<-window.names
names(site.counts)<-c("scaffold","pos","AA","AB","BB","NN")

#Write and save the table
write.table(site.counts, file="8C_site.counts_hetdev=0.2+WS50K.txt",sep="\t",quote=FALSE,row.names=FALSE)
