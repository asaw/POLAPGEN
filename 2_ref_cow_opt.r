#The best waping for reference chromatograms selection

mv=c(132,176,220,275)
tv=c(10,12,15,18)

track_data_general=paste("/home/users/asawigr/reef/zad15/Dane/3_wysiew_liscie/",sep="")
track_results=paste("/home/users/asawigr/reef/zad15/3_wysiew_liscie/warping",sep="")
setwd(track_data_general

#Table of observation import
LAB=read.csv("tabela_obserwacji_zad15_doœw3.csv")
obk=unique(LAB[,2])
ilm=length(obk)
setwd(track_results)
ind_ref=read.table("id_line_ref.txt")
ind_ref=as.real(ind_ref)
ind_ref
print(obk[ind_ref]) #reference line/variety
REF=as.matrix(read.table(paste("REF_",id_wave,".txt",sep="")))
ref=REF[ind_ref,]
time=as.matrix(read.table("time.txt"))
L=13201
ptl=Sys.time()
print(ptl)
time=time[1:L]


###########################
###function produces warping effect

fun_wef=function(L,time,ref,REF,ind_ref,mv,tv,id_wave,ilm){
time=time[1:L]
filenames <- dir(pattern=paste("gridref_",as.character(id_wave),"_",sep=""),full.names=F)  
P <- lapply(filenames, read.table,header=T,fill=T)
library(Matrix)
ind_ref #reference chromatogram number
###########################
T=ref
Pdiff=list(1:(ilm-1))
for (dd in 1:ilm){
if (dd!=ind_ref){
Pdiff[[dd]]=REF[dd,1:L]
}
}
############################
#Warping effect

Lp=(L-1)
R=NULL
count=0
bestWarp=0
for (countm in 1:length(mv)){
for (countt in 1:length(tv)){ 
m=mv[countm]
t=tv[countt]
X=NULL
print(m)
print(t)
if (m>(t+3)) {
count=count+1
X=T
for (ee in 1:ilm){
if (ee!=ind_ref){
X=rbind(X,P[[ee]][count,3:(L+2)])
}
}
Y=X

X=X/norm(t(t(X)),"F")
S=svd(X)
Simp=S$d[1]^4+S$d[2]^4

n=list(1:(2*(ilm-1)))

for (ff in 1:ilm){
if (ff!=ind_ref){
n[[(2*ff-1)]]=P[[ff]][count,3:(L+2)]
n[[(2*ff-1)]]=t(t(n[[(2*ff-1)]]))
n[[(2*ff-1)]]=norm(n[[(2*ff-1)]],"F")
n[[2*ff]]=Pdiff[[ff]]
n[[2*ff]]=t(t(n[[2*ff]]))
n[[2*ff]]=norm(n[[2*ff]],"F")
}
}

c=list(1:(ilm-1))

for (ff in 1:ilm){
if (ff!=ind_ref){
c[[ff]]=abs((n[[(2*ff-1)]]-n[[2*ff]])/n[[2*ff]])
}
}

Peak=0

for (ff in 1:ilm){
if (ff!=ind_ref){
Peak=Peak+(1-(min(c[[ff]],1))^2)
}
}

Peak=Peak/(ilm-1)

Warp=Simp+Peak


row=cbind(m,t,Simp,Peak,Warp)
print(row)
R=rbind(R,row)

if (bestWarp<Warp){
bestWarp=Warp
bestm=m
bestt=t
PCAmatrix=Y
}

}
}
}
print(bestWarp)
print(bestm)
print(bestt)


write.table(R,paste("warp_eff_grid_ref_",as.character(id_wave),".txt",sep=""))
write.table(PCAmatrix,paste("PCAmatrix_grid_ref_",as.character(id_wave),".txt",sep="")) #PCAmatrix - data after the best warping (without normalization as in X)


#PCAmatrix, X, Y has order as ref is first, others: (P, P2, Pdiff) has order as labels
P2=REF

pdf2=paste("ref_",id_wave,"_before_cow.pdf",sep="")
pdf(pdf2,width=10,height=6)
title2=paste("reference chromatograms after differentiation, before COW - ",id_wave," nm",sep="")
plot(time,P2[1,],type="l",ylim=c(-0.06,0.06),main=title2, xlab="Time [min]",ylab="Differented absorbance",col=1)
for (cc in 1:ilm){
lines(time,P2[cc,],col=cc)
}
dev.off()

# PCAnormal - matrix with data from PCAmatrix, but with order as in labels
if (ind_ref==1) {PCAnormal=PCAmatrix} else {if (ind_ref==ilm) {PCAnormal=rbind(PCAmatrix[2:ilm,],PCAmatrix[1,])} else {PCAnormal=rbind(PCAmatrix[2:ind_ref,],PCAmatrix[1,],PCAmatrix[(ind_ref+1):ilm,])}}

pdf3=paste("ref_",id_wave,"_after_cow.pdf",sep="")
title3=paste("reference chromatograms after differentiation, after COW - ",id_wave," nm",sep="")
pdf(pdf3,width=10,height=6)
plot(time,PCAnormal[1,],ylim=c(-0.06,0.06),type="l",main=title3, xlab="Time [min]",ylab="Differented absorbance",col=1)
for (cc in 1:ilm){
lines(time,PCAnormal[cc,],col=cc)
}
dev.off()


#########standard

pdf2=paste("ref_",id_wave,"_before_cow_standard.pdf",sep="")
pdf(pdf2,width=10,height=6)
title2=paste("reference chromatograms after differentiation, before COW - ",id_wave," nm",sep="")
plot(time,P2[1,],type="l",ylim=c(-0.3,0.3),xlim=c(12,12.5),main=title2, xlab="Time [min]",ylab="Differented absorbance",col=1)
for (cc in 1:ilm){
lines(time,P2[cc,],col=cc)
}
dev.off()

# PCAnormal - matrix with data from PCAmatrix, but with order as in labels
if (ind_ref==1) {PCAnormal=PCAmatrix}
if (ind_ref==ilm) {PCAnormal=rbind(PCAmatrix[2:ilm,],PCAmatrix[1,])}
PCAnormal=rbind(PCAmatrix[2:ind_ref,],PCAmatrix[1,],PCAmatrix[(ind_ref+1):ilm,])

pdf3=paste("ref_",id_wave,"_after_cow_standard.pdf",sep="")
title3=paste("reference chromatograms after differentiation, after COW - ",id_wave," nm",sep="")
pdf(pdf3,width=10,height=6)
plot(time,PCAnormal[1,],ylim=c(-0.3,0.3),xlim=c(12,12.5),type="l",main=title3, xlab="Time [min]",ylab="Differented absorbance",col=1)
for (cc in 1:ilm){
lines(time,PCAnormal[cc,],col=cc)
}
dev.off()


print("The correct ending")
}

id_wave=as.real(id_wave)
fun_wef(L,time,ref,REF,ind_ref,mv,tv,id_wave,ilm)

