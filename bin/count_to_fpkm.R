expr1<-read.delim('count1.txt',header = T)
expr1<-expr1[!duplicated(expr1$Geneid),]
mygeneLeng<-expr1$Length
names(mygeneLeng)<-expr1$Geneid

mycount<-expr1[,7:ncol(expr1)]
rownames(mycount)<-expr1$Geneid
mycount<-mycount

mycount<-sweep(mycount,MARGIN = 2,colSums(mycount),'/')

myfpkm<-sweep(mycount,MARGIN = 1,mygeneLeng,'/')*10^9

write.table(cbind(id=rownames(myfpkm),myfpkm),'fpkm_result.txt',row.names = F,sep = '\t',quote = F)
