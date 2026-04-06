#FGA calculation####
setwd("./job_example/")

# import CNV data for each organoid####
file_names <- list.files(pattern="cns.txt$")
data_merged <- read.table(file_names[1],header=T,sep = "\t")
data_merged$ID=file_names[1]
table(data_merged$Otherinfo)
data_merged=subset(data_merged,data_merged$Otherinfo!="")
data_merged$end=as.numeric(data_merged$end)
data_merged$start=as.numeric(data_merged$start)
data_merged$length=data_merged$end-data_merged$start
FGA=aggregate(data_merged$length,by=list(type=data_merged$ID),sum)
FGA$ratio=FGA$x/3000000000
data_merged <- read.table(file_names[1],sep="\t",header=T)
data_merged$ID=file_names[1]

for (i in 2:length(file_names)) {
  data_temp <- read.table(file_names[i],sep="\t",header=T)
  data_temp$ID=file_names[i]
  data_merged <- rbind(data_merged, data_temp)
}
data_merged=data_merged[,c(8,1:7)]
data_merged$ID=gsub("[.].*$","",data_merged$ID)

data_merged$length=data_merged$end-data_merged$start
FGA=aggregate(data_merged$length,by=list(type=data_merged$ID),sum)
FGA$ratio=FGA$x/3000000000
FGA$type=gsub("[_].*$","",FGA$type)
write.csv(FGA, "FGA_organoid.csv", row.names = FALSE)