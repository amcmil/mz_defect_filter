#m/z defect filtering commands
#12-04-2016

#--------------------------------------------------------------------------------
#The following chunck of code will remove features in your LC-MS data identified as salt clusters based on m/z defect, retention time and an inclusion list derived from the human metabolome database (included in this package).
#your feature table must have LC-MS features in rows, m/z in column 1 (labelled "mz"), followed by retention time in column 2 (labelled "rt").
#-------------------------------------------------------------------------------- 

#read in feature table you want to filter
mydata<-read.table("example_data_pos.txt",  header=T, check.names=F, row.names=1, sep="\t")

#set rt threshold
rt=0.8

#set mass accuracy (in Daltons) for inclusion list matching
ma=0.01

#calculate m/z defect using modulus function
out<-matrix(data=NA, nrow =nrow(mydata), ncol=1)
rownames(out)<-rownames(mydata)

for(i in 1:nrow(mydata)){
out[i,]<-mydata[i,1]%%1
colnames(out)<-"md" 
}

mydata<-cbind(mydata,out)

#fit mydata to m/z defect maxima equation
f<- function(x){
y<-0.00112*x + 0.01953
return(y)
}

fit<-as.data.frame(t(apply(mydata,1,f)))

#read in hmdb inclusion list
hmdb<-read.table("hmdb_inclusions_list_pos.txt",  header=T, check.names=F, row.names=1, sep="\t")

#find masses in inclusion list in mydata 
in_list<-hmdb[which(hmdb$"inclusion"=="y"),]
cmp <- function(mydata, in_list, cutoff=ma){abs(in_list-mydata) <= cutoff}
match<-which(outer(in_list$mz, mydata$mz, cmp), arr.ind=TRUE)
fc <- factor(match[,2])
I<-rownames(mydata)[as.numeric(levels(fc))]

#keep rows in mydata if m/z defect is less than or equal to fitted value OR rt > threshold OR value is in inclusion list

filtered<-mydata[which((mydata$"md"<= fit$"mz")|rownames(mydata) %in% I | mydata$"rt" >rt),]

#write filtered table to file
write.table(filtered,"feature_list_mz_defect_filtered.txt",sep="\t",col.names=NA)

#plot mz vs md of raw and filtered features 
pdf("mz defect plot before and after filter.pdf", height=6.5,width=5)  
par(mfrow=c(2,1),mai=c(0.9,0.9,0.4,0.5))
plot(mydata$"mz",mydata$"md",cex.axis=0.8,col=adjustcolor("black",alpha=0.3),pch=20,cex=0.8,ylim=c(0,1),xlim=c(50,750),ylab="md",xlab="m/z",main="Raw",cex.lab=0.8,cex.main=0.8)
plot(filtered$"mz",filtered$"md",cex.axis=0.8,col=adjustcolor("black",alpha=0.3),pch=20,cex=0.8,ylim=c(0,1),xlim=c(50,750),ylab="md",xlab="m/z",main="Filtered",cex.lab=0.8,cex.main=0.8)
dev.off()
