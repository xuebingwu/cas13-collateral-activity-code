## fig1h ###

# load fig1h data 
x=read.csv("./fig1h.csv",row.name=1)

# filter out genes never expressed in any sample
x <- x[apply(x,1,max) > 0,]

# specific genes
DsRed=grep('DsRed',rownames(x))
Cas13=grep('RfxCas13d',rownames(x))
BFP=grep('TagBFP',rownames(x))
mito=grep('^MT-',rownames(x))
mouse=grep('mouse',rownames(x))
gapdh=which(rownames(x) == 'GAPDH')
actb=which(rownames(x) == 'ACTB')
hn=which(rownames(x) == 'HNRNPA2B1')
cd=which(rownames(x) == 'CD99')
b4=which(rownames(x) == 'B4GALNT1')

# normalize by mouse spike in

for (i in 1:(ncol(x)-1)){
  x[,i] = x[,i]/ sum(x[mouse,i])*1e6 
}

# average over replicates
KDg3  <- apply(x[,1:2],1,mean) # DsRed KD
NTg11 <- apply(x[,3:4],1,mean) # Control


pdf("fig1h.pdf",width=4,height=4)
par(pty="s")

ctrl=NTg11
knockdown=KDg3

plot(ctrl,knockdown,cex=0.1,pch=16,xlim=c(1,max(ctrl)),ylim=c(1,max(knockdown)),log='xy',bty='l',col=rgb(0,0,0,0.1),xlab='Expression in TPM (Control)\n',ylab='\nExpression in TPM (Knockdown)',main='')

points(ctrl[mito],knockdown[mito],pch=16,col='cyan',cex=0.5)
points(ctrl[mouse],knockdown[mouse],pch=16,col=rgb(1,0,1,0.2),cex=0.1)
points(ctrl[DsRed],knockdown[DsRed],pch=16,col='red')
points(ctrl[Cas13],knockdown[Cas13],pch=16,col='green')
points(ctrl[BFP],knockdown[BFP],pch=16,col='blue')
points(ctrl[gapdh],knockdown[gapdh],pch=16,col='orange')
points(ctrl[actb],knockdown[actb],pch=16,col='purple')

dev.off()


