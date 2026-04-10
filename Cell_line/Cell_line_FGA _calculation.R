###Cell line_FGA###
CNA=read.csv("./job_example/CCLE_wes_segment_cn.csv",sep=",")
info=read.csv("./job_example/sample_info.csv")
CNA$size=CNA$End-CNA$Start
CNA=CNA[which(abs(CNA$Segment_Mean) > 0.2 & CNA$Status !=0),]
b=aggregate(CNA$size, by=list(type=CNA$DepMap_ID),sum)
info=info[,c(1,3,13)]
colnames(info)[2]="name"
colnames(info)[3]="type"
table(info$type)
info=subset(info,info$type=="Colon/Colorectal Cancer")
colnames(info)="ID"
colnames(b)="ID"
mer=merge(info,b,by="ID",all.x=T)
colnames(mer)[2]="name"
colnames(mer)[4]="value"
mer=mer[,c(2,4)]
mer$ratio=mer$value/3000000000
mer=subset(mer,!is.na(mer$ratio))
write.csv(mer,"./job_example/FGA.csv")
