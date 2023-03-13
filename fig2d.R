
# load DsRed data
x0=read.csv("./fig1g.csv",row.name=1)
# filter out genes never expressed in any sample
x0 <- x0[apply(x0,1,max) > 0,]
# assign genes to classes
x0$class='Human' # default
x0$class[grep('plasmid',rownames(x0))] = 'Plasmid'
x0$class[grep('DsRed',rownames(x0))] = 'DsRed'
x0$class[grep('RfxCas13d',rownames(x0))]='Cas13'
x0$class[grep('TagBFP',rownames(x0))]='BFP'
x0$class[grep('^MT-',rownames(x0))]='Mito'
x0$class[grep('mouse',rownames(x0))]='Mouse'
for (i in 1:(ncol(x0)-1)){
  x0[,i] = x0[,i]/ sum(x0[x0$class=='Mito',i])*1e6 
}
# average over replicates
KDg3  <- apply(x0[,2:3],1,mean) # DsRed KD
NTg11 <- apply(x0[,5:6],1,mean) # Control

# load endogenous gene data
x=read.csv("./fig2d.csv",row.name=1)

# filter out genes never expressed in any sample
x <- x[apply(x,1,max) > 0,]

# assign genes to classes
x$class='Human' # default
x$class[grep('^MT-',rownames(x))]='Mito'

# normalize to mito rnas
for (i in 1:(ncol(x)-1)){
  x[,i] = x[,i]/ sum(x[x$class=='Mito',i])*1e6 
}

# average over replicates
NT <- apply(x[,1:2],1,mean)  
HNRNPA2B1 <- apply(x[,3:4],1,mean)  
CD99 <- apply(x[,5:6],1,mean)  
B4GALNT1 <- apply(x[,7:8],1,mean)  

pdf("fig2d.pdf",width=4,height=4)

par(mfrow=c(5,1),mar=c(0,1,0,1),bty='l',cex=1)
xlim=c(-1.5,1.5)
ylim=c(0,2)

sub= x0$class=='Human' & NTg11 > 0
mito = x0$class=='Mito' & NTg11 > 0
plot(0,xlim=xlim,ylim=ylim,axes=FALSE, ann=FALSE, frame.plot=FALSE,cex=0,col='blue')
polygon(density(log2(KDg3[sub]/NTg11[sub]),n=2048), col=rgb(0,0.4,1,0.5),border=NA)
lines(density(log2(KDg3[mito]/NTg11[mito])),col='gray')

med= median(log2(KDg3[sub]/NTg11[sub])) - median(log2(KDg3[mito]/NTg11[mito]))
d=density(log2(KDg3[mito]/NTg11[mito]))
ks <- ks.test(log2(KDg3[sub]/NTg11[sub]),log2(KDg3[mito]/NTg11[mito]))
legend(-1.8,2.2,legend=paste(round(2^med-1,digits=2)*100,'%, P=',signif(ks$p.value,digits=1),sep=''),bty='n')
legend('right',legend='DsRed',bty='n')
abline(v=0,lty=2)

text(-0.8,1.0,'Other RNA',col=rgb(0,0.4,1,0.5))
text(0.05,1.4,'mt-RNA',col='gray')

sub= x$class=='Human' & NT > 0
mito= x$class=='Mito' & NT > 0
plot(0,xlim=xlim,ylim=ylim,axes=FALSE, ann=FALSE, frame.plot=FALSE,cex=0,col='blue')
polygon(density(log2(HNRNPA2B1[sub]/NT[sub]),n=2048), col=rgb(0,0.4,1,0.5),border=NA)
lines(density(log2(HNRNPA2B1[mito]/NT[mito])),col='gray')
med=median(log2(HNRNPA2B1[sub]/NT[sub])) - median(log2(HNRNPA2B1[mito]/NT[mito])) 
ks <- ks.test(log2(HNRNPA2B1[sub]/NT[sub]),log2(HNRNPA2B1[mito]/NT[mito]))
legend(-1.8,2.2,,legend=paste(round(2^med-1,digits=2)*100,'%, P=',signif(ks$p.value,digits=1),sep=''),bty='n')
legend('right',legend='HNRNPA2B1',bty='n')
abline(v=0,lty=2)


plot(0,xlim=xlim,ylim=ylim,axes=FALSE, ann=FALSE, frame.plot=FALSE,cex=0,col='blue')
polygon(density(log2(CD99[sub]/NT[sub]),n=2048), col=rgb(0,0.4,1,0.5),border=NA)
lines(density(log2(CD99[mito]/NT[mito])),col='gray')
med=median(log2(CD99[sub]/NT[sub])) - median(log2(CD99[mito]/NT[mito])) 
ks <- ks.test(log2(CD99[sub]/NT[sub]),log2(CD99[mito]/NT[mito]))
legend(-1.8,2.2,,legend=paste(round(2^med-1,digits=2)*100,'%, P=',signif(ks$p.value,digits=1),sep=''),bty='n')
legend('right',legend='CD99',bty='n')
abline(v=0,lty=2)


plot(0,xlim=xlim,ylim=ylim,axes=FALSE, ann=FALSE, frame.plot=FALSE,cex=0,col='blue')
polygon(density(log2(B4GALNT1[sub]/NT[sub]),n=2048), col=rgb(0,0.4,1,0.5),border=NA)
lines(density(log2(B4GALNT1[mito]/NT[mito])),col='gray')
med=median(log2(B4GALNT1[sub]/NT[sub])) - median(log2(B4GALNT1[mito]/NT[mito]))
ks <- ks.test(log2(B4GALNT1[sub]/NT[sub]),log2(B4GALNT1[mito]/NT[mito]))
legend(-1.8,2.2,,legend=paste(round(2^med-1,digits=2)*100,'%, P=',signif(ks$p.value,digits=1),sep=''),bty='n')
legend('right',legend='B4GALNT1',bty='n')
abline(v=0,lty=2)

axis(1,at= c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5),c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5))

plot(0,xlim=xlim,ylim=ylim,axes=FALSE, ann=FALSE, frame.plot=FALSE,cex=0,col='blue')
text(0,0.6,'Log2 fold change')

dev.off()

