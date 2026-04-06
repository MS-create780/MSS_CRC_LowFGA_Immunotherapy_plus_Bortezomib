#FGA####
setwd("./job_example/")
####import CNV data for each patient###
filenames=dir()
data=data.frame()
for (i in filenames){
  path=paste0(getwd(),"\\",i)
  data=rbind(data,read.table(path,header=T))
}
CNA=data
CNA$mer=paste(CNA$chrom,CNA$start,CNA$end)
id=as.data.frame(table(CNA$sample_id))
list=id$Var1
r=as.data.frame(c(1:16))
r=as.data.frame(t(r))
colnames(r)=colnames(CNA)
for(i in list){
  p=subset(CNA,CNA$sample_id=="WJ_1")
  p2=p[-which(duplicated(p$mer)),]
  r=rbind(r,p2)
}
z=as.data.frame(table(r$sample_id))
r=r[-1,]
r=r[r$copy_number!=2,]
CNA=r
CNA$size=CNA$end-CNA$start
b=aggregate(CNA$size, by=list(type=CNA$sample_id),sum)
b$ratio=b$x/3000000000
write.csv(b,"./job_example/CRC_sample_FGA.csv")