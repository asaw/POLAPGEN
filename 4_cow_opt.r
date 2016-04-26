#The best waping for variety selection

####################The number of considered variety identifing on the list
#id_variety="072" #change od 1 do 105(length(obk))
#id_wave="330" #change
#variety_ref=42
mv=c(132,176,220,275)
tv=c(10,12,15,18)

track_data_wave=paste("/home/users/asawigr/reef/zad15/Dane/3_wysiew_liscie/",id_wave,"nm",sep="")
track_data_general=paste("/home/users/asawigr/reef/zad15/Dane/3_wysiew_liscie/",sep="")
track_results=paste("/home/users/asawigr/reef/zad15/3_wysiew_liscie/warping",sep="")

setwd(track_results)
variety_ref=read.table("id_line_ref.txt")  #the reference variety number loading
variety_ref=as.real(variety_ref)
variety_ref

REF_info=read.table("REF_info.txt") #the reference chromatogram number in each variety loading
ind_ref=as.real(REF_info[as.real(id_variety),1])
ind_ref

setwd(track_data_wave) 
filenames <- dir(pattern="arw",full.names=F)  
matrices <- lapply(filenames, read.table,header=T,fill=T)
ilm=length(matrices) 
time=as.vector(matrices[[1]][2402:15603,1])
L=13201
ptl=Sys.time()
print(ptl)
ilm


#########################
setwd(track_data_general)
#Table of observation import
LAB=read.csv("tabela_obserwacji_zad15_doœw3.csv")

obk=unique(LAB[,2])
#write.csv(obk,"lista_obiektow.csv")

#numkol - numbers from the data table set in order of loading files
#samnr - variety name according to the order in the data table
#nrtab - numbers successively loaded files on the selected variety
nrtab=NULL
Labels=NULL
countm=0
samnr=as.character(LAB[,2])
for (ll1 in 1:ilm){
numkol=as.real(substr(as.character(matrices[[ll1]][1,2]),1,4))
if (samnr[numkol]==obk[as.real(id_variety)]){
nrtab=rbind(nrtab,ll1)
#Sample labels creation - reference chromatogram is not as a first, it is just on its place
Labels=rbind(Labels,LAB[numkol,])
#Samples from selected variety counting
countm=countm+1
}
}

print(countm)
ilm=countm

#matrices will include only selected variety
setwd(track_data_wave) 
matrices <- lapply(filenames[nrtab], read.table,header=T,fill=T)
setwd(track_results) 
#text00=paste("labels_samp_",id_variety,"_",id_wave,".txt",sep="")

#write.table(Labels,text00)


########Data normalization
#All chromatograms are dividing by the mass of each chromatogram
masa=NULL
for (mi in 1:ilm){
masa=c(masa,matrices[[mi]][1,3])
}

for (mi in 1:ilm){
matrices[[mi]][2402:15603,2]=matrices[[mi]][2402:15603,2]/masa[mi]
}

###Rys.a.
pdf1=paste(id_variety,"_",id_wave,"_after_normalization_check.pdf",sep="")
title1=paste(id_variety," after normalization - ",id_wave," nm",sep="")
pdf(pdf1,width=10,height=6)
plot(time,as.vector(matrices[[1]][2402:15603,2]),type="n",ylim=c(0,3),main=title1,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices[[cc]][2402:15603,2]),col=cc)
}
dev.off()

########Baseline removal based on differentiation:
for (mi in 1:ilm){
matrices[[mi]][2402:15602,2]=diff(matrices[[mi]][2402:15603,2])
}

###Rys.b.

pdf2a=paste(id_variety,"_",id_wave,"_before_cow_check.pdf",sep="")
title2a=paste(id_variety," before_cow_check - ",id_wave," nm",sep="")
pdf(pdf2a,width=10,height=6)
plot(time,as.vector(matrices[[1]][2402:15603,2]),type="n",ylim=c(-0.06,0.06),main=title2a,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices[[cc]][2402:15603,2]),col=cc)
}
dev.off()




pdf2a=paste(id_variety,"_",id_wave,"_before_cow_standard_check.pdf",sep="")
title2a=paste(id_variety," befpre_cow_check - ",id_wave," nm",sep="")
pdf(pdf2a,width=10,height=6)
plot(time,as.vector(matrices[[1]][2402:15603,2]),type="n",ylim=c(-0.3,0.3),xlim=c(12,12.5),main=title2a,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices[[cc]][2402:15603,2]),col=cc)
}
dev.off()



setwd(track_results)
time=time[1:L]

#########################
#ref importing (after all reference chromatograms warping amoung all varieties)
REF1=read.table(paste("PCAmatrix_grid_ref_",id_wave,".txt",sep=""))

#if variety_ref is considered reference variety (variety_ref - the reference variety number taken from the order obk), then its profile after moving is at the first place,
#otherwise (if another variety) indexes have to be moved

if (as.real(id_variety)==variety_ref){move_id=1} else {
if (as.real(id_variety)<variety_ref){move_id=as.real(id_variety)+1}
if (as.real(id_variety)>variety_ref){move_id=as.real(id_variety)}
}
ref=as.real(REF1[(move_id),])
#########################


#################################################
time=time[1:L]
pat1=paste("grid",id_variety,"_",id_wave,"_",sep="")
filenames <- dir(pattern=pat1,full.names=F)  
P <- lapply(filenames, read.table,header=T,fill=T)
library(Matrix)
ind_ref #the reference chromatogram number 
###########################

T=ref
Pdiff=list(1:(ilm-1))
for (dd in 1:ilm){
if (dd!=ind_ref){
Pdiff[[dd]]=as.vector(matrices[[dd]][2402:15602,2])
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

table1=paste("warp_eff_grid_",id_wave,"_",id_variety,".txt",sep="")
write.table(R,table1)
table2=paste("PCAmatrix_grid_",id_wave,"_",id_variety,".txt",sep="")
write.table(PCAmatrix,table2) #PCAmatrix - data after the best warping (without normalization as in X)



############Rys.b. i c.
T=ref
P2=NULL
Pdiff=list(1:(ilm-1))
for (dd in 1:ilm){
P2=rbind(P2,as.vector(matrices[[dd]][2402:15602,2]))
if (dd!=ind_ref){
Pdiff[[dd]]=as.vector(matrices[[dd]][2402:15602,2])
}
}
P2[ind_ref,]=ref
############################
#table2=paste("PCAmatrix_grid_",id_variety,"_",id_wave,".txt",sep="")
#PCAmatrix=read.table(table2) #PCAmatrix - data after the best warping (without normalization as in X)



#PCAmatrix, X, Y has order as ref is first, others: (P, P2, Pdiff) has order as labels

pdf2=paste(id_variety,"_",id_wave,"_before_cow.pdf",sep="")
pdf(pdf2,width=10,height=6)
title2=paste(id_variety," after differentiation, before COW - ",id_wave," nm",sep="")
plot(time,P2[1,],type="l",ylim=c(-0.06,0.06),main=title2, xlab="Time [min]",ylab="Differented absorbance",col=1)
for (cc in 1:ilm){
lines(time,P2[cc,],col=cc)
}
dev.off()

# PCAnormal - matrix with data from PCAmatrix, but with order as in labels
if (ind_ref==1) {PCAnormal=PCAmatrix}
if (ind_ref==ilm) {PCAnormal=rbind(PCAmatrix[2:ilm,],PCAmatrix[1,])}
PCAnormal=rbind(PCAmatrix[2:ind_ref,],PCAmatrix[1,],PCAmatrix[(ind_ref+1):ilm,])

pdf3=paste(id_variety,"_",id_wave,"_after_cow.pdf",sep="")
title3=paste(id_variety," after differentiation, after COW - ",id_wave," nm",sep="")
pdf(pdf3,width=10,height=6)
plot(time,PCAnormal[1,],ylim=c(-0.06,0.06),type="l",main=title3, xlab="Time [min]",ylab="Differented absorbance",col=1)
for (cc in 1:ilm){
lines(time,PCAnormal[cc,],col=cc)
}
dev.off()


#########standard

pdf2=paste(id_variety,"_",id_wave,"_before_cow_standard.pdf",sep="")
pdf(pdf2,width=10,height=6)
title2=paste(id_variety," after differentiation, before COW - ",id_wave," nm",sep="")
plot(time,P2[1,],type="l",ylim=c(-0.3,0.3),xlim=c(12,12.5),main=title2, xlab="Time [min]",ylab="Differented absorbance",col=1)
for (cc in 1:ilm){
lines(time,P2[cc,],col=cc)
}
dev.off()

# PCAnormal - matrix with data from PCAmatrix, but with order as in labels
if (ind_ref==1) {PCAnormal=PCAmatrix}
if (ind_ref==ilm) {PCAnormal=rbind(PCAmatrix[2:ilm,],PCAmatrix[1,])}
PCAnormal=rbind(PCAmatrix[2:ind_ref,],PCAmatrix[1,],PCAmatrix[(ind_ref+1):ilm,])

pdf3=paste(id_variety,"_",id_wave,"_after_cow_standard.pdf",sep="")
title3=paste(id_variety," after differentiation, after COW - ",id_wave," nm",sep="")
pdf(pdf3,width=10,height=6)
plot(time,PCAnormal[1,],ylim=c(-0.3,0.3),xlim=c(12,12.5),type="l",main=title3, xlab="Time [min]",ylab="Differented absorbance",col=1)
for (cc in 1:ilm){
lines(time,PCAnormal[cc,],col=cc)
}
dev.off()

#Pdiff_ord has the order of rows as in Diff1, PCAmatrix, X, Y,
#so as Labels, that ref is at the first place
#Pdiff_ord=rbind(P2[ind_ref,],P2[1:(ind_ref-1),],P2[(ind_ref+1):ilm,])

#pdf2=paste(id_variety,"_",id_wave,"_before_cow.pdf",sep="")
#pdf(pdf2,width=10,height=6)
#title2=paste(id_variety," after differentiation, before COW - ",id_wave," nm",sep="")
#plot(time,Pdiff_ord[1,],type="l",ylim=c(-0.3,0.3),main=title2, xlab="Time [min]",ylab="Differented absorbance",col=1)
#for (cc in 1:ilm){
#lines(time,Pdiff_ord[cc,],col=(cc+1))
#}
#dev.off()


#pdf3=paste(id_variety,"_",id_wave,"_po_cow.pdf",sep="")
#title3=paste(id_variety," after differentiation, after COW - ",id_wave," nm",sep="")
#pdf(pdf3,width=10,height=6)
#plot(time,PCAmatrix[1,],ylim=c(-0.3,0.3),type="l",main=title3, xlab="Time [min]",ylab="Differented absorbance",col=1)
#for (cc in 1:ilm){
#lines(time,PCAmatrix[cc,],col=(cc+1))
#}
#dev.off()

print("The correct ending")

