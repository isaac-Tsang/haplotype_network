library(pacman) 
p_load(pegas,bruceR,"ggsci",scales,"gridExtra")
set.wd("H:/paper/GLR/约稿综述/frontier投稿/其他老师修改意见/黎志康老师/predominant gcHap")
#rm(list = ls())
#dev.off()
Pop <- import('3k-3024分群.xlsx',as="tibble",header = F)
colnames(Pop) = c("SNPID", "Subgroup") 

Pop <-  Pop %>% filter(Subgroup != "Adm" )

## gene =========================
gene=c(
  "Os02t0787600-01.CDS"
)


## prepare----------------------------------------------------------------------

myPrepare <- function(i,Num){
  aa=gene[i]
  file1=as.character(paste(aa,"/",aa,".geno", sep=""))
  df <- fread(file1,header = F)
  df1 <- as.data.frame(t(df[,-c(1:7)]))
  
  #df1[,2] <- gsub("-", "T", df1[,2])
 
  df2 <- df1 %>% unite(Seq, !V1, sep = "") 
  
  df2 <- df2 %>% mutate(Seq = gsub("-", "T", Seq))
 
  colnames(df2) <- c("SNPID",  "seq")
 
  y <- df2 %>% inner_join(Pop,by="SNPID") %>% filter(Subgroup!="na")
  
  head(y)
  
  a=as.data.frame(table(y$seq))
  tail(arrange(a,Freq),n=20)
  b=as.data.frame(tail(arrange(a,Freq),n=20))
  c=b %>% filter(Freq >Num)
  c
  colnames(c)=c("seq","Freq")
  c
  d = merge(y,c,by = "seq")
  head(d)
 
  e=d[,c(2,3,1)]
  
  
  head(e)
  
  file2=as.character(paste(aa,"/",aa,"-",Num,".CSV", sep=""))
  write.table(e,file2,sep = ",",quote = F,row.names = F,col.names = F)
  
  f = d[, c("SNPID", "seq")]
  fasta_seq <- paste("> ", e$SNPID, "\n", e$seq, "\n", sep = "")
  file11 = paste(aa, "/", aa, "-", Num, ".fasta", sep = "")
  writeLines(fasta_seq, file11)
  
  rm(a,b,c,d,e,y,fasta_seq)
}
Fre=50
S=1
myPrepare(i=S,Num=Fre)
aa=gene[S]

## start plot-------------------------------------------------------------------

file3=as.character(paste(aa,"/",aa,"-",Fre,".CSV", sep=""))
y <- read.table(file3, sep = ",", header = F)
names(y) <- c("SNPID", "Subgroup", "seq")
head(y)
a=as.data.frame(table(y$seq))
arrange(a,Freq)

file4=as.character(paste(aa,"/",aa,"-",Fre,".fasta", sep=""))

x <- read.dna(file4, format = "fasta")

head(x)

h <- haplotype(x)

print(h, T)


#输出对应Hap的文件，显示有几个就输出几个
qq=attributes(h)

print(qq)

qa=qq$index[[1]]
y1=y[qa,]
write.csv(y1,"Hap1.csv",row.names=FALSE)
qb=qq$index[[2]]
y2=y[qb,]
write.csv(y2,"Hap2.csv",row.names=FALSE)
qc=qq$index[[3]]
y3=y[qc,]
write.csv(y3,"Hap3.csv",row.names=FALSE)
qd=qq$index[[4]]
y4=y[qd,]
write.csv(y4,"Hap4.csv",row.names=FALSE)
#qe=qq$index[[5]]
#y5=y[qe,]
#write.csv(y5,"Hap5.csv",row.names=FALSE)
#qb=qq$index[[6]]
#y6=y[qb,]
#write.csv(y6,"Hap6.csv",row.names=FALSE)
#qc=qq$index[[7]]
#y7=y[qc,]
#write.csv(y7,"Hap7.csv",row.names=FALSE)
#qd=qq$index[[8]]
#y8=y[qd,]
#write.csv(y8,"Hap8.csv",row.names=FALSE)
#qe=qq$index[[9]]
#y9=y[qe,]
#write.csv(y9,"Hap9.csv",row.names=FALSE)





attributes(h)
net <- haploNet(h)
snp_type <- stack(setNames(attr(h, "index"), rownames(h)))
result <- cbind(snp_type, y[snp_type[[1]],])

ind.hap2<-with(
  stack(setNames(attr(h, "index"), rownames(h))), 
  table(hap=ind, pop=y[[2]][values])
)

par(bg="white")
par(mar=c(1,1,4,1))

#c1=c("#00bdcd","#fadc82","#006a7b","#ef1828","#f88421")#HZQ
#c1 = c("#4daf4a", "#ff7f00", "#984ea3", "#FFA07A") #isaac
#c1 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442") #isaac
c1 <- c("#DAB6C4", "#957DAD", "#27AE60", "#3498DB")

show_col(c1)

plot(net, size=attr(net, "freq"), scale.ratio = 240, asp = 1,
     pie=ind.hap2,
     bg=c1, 
     legend = FALSE,  fast = FALSE, show.mutation = 1, 
     threshold = c(2,3,3), 
     labels = TRUE,
     col.link = "black", cex = 1, lwd = 1, lty = 1, font=2)

aa
print.default(net)
ind.hap2


legend(-600,
       1400,
       colnames(ind.hap2), 
       #col=col.vec, 
       col=c1, 
       pch=15, 
       cex = 1.2, border ="white", box.col = "white",
       y.intersp=0.8,x.intersp = 0.5)
title(main = "Haplotype Network Plot of Os02t0787600-01.CDS",
      line = -0.1, col.main= "black", cex.main=0.8)

qq=attributes(h)
qa=qq$index[[1]]
y1=y[qa,]
#View(y1)
table(y1$seq)
qa=qq$index[[2]]
y1=y[qa,]
table(y1$seq)
qa=qq$index[[3]]
y1=y[qa,]
table(y1$seq)
qa=qq$index[[4]]
y1=y[qa,]
table(y1$seq)
qa=qq$index[[5]]
y1=y[qa,]
table(y1$seq)
qa=qq$index[[6]]
y1=y[qa,]
table(y1$seq)
qa=qq$index[[7]]
y1=y[qa,]
table(y1$seq)
qa=qq$index[[8]]
y1=y[qa,]
table(y1$seq)
#qa=qq$index[[9]]
#y1=y[qa,]
#table(y1$seq)
#qa=qq$index[[10]]
#y1=y[qa,]
#table(y1$seq)
