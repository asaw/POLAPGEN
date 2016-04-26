#reference chromatogram plots after cow and back to not differenced data
id_wave_ref=280
id_ch=2

mv=c(132)
tv=c(18)

track_data_general=paste("/home/users/asawigr/reef/zad15/Dane/3_wysiew_liscie/",sep="")
track_results=paste("/home/users/asawigr/reef/zad15/3_wysiew_liscie/warping",sep="")

#######################
setwd(track_data_general)
#Table of observation import
LAB=read.csv("tabela_obserwacji_zad15_doœw3.csv")
obk=unique(LAB[,2])


#############
id_wave=280
id_wave=as.character(id_wave)
track_data_wave=paste("/home/users/asawigr/reef/zad15/Dane/3_wysiew_liscie/",id_wave,"nm",sep="")
setwd(track_data_wave) 
filenames1 <- dir(pattern="arw",full.names=F)  
matrices <- lapply(filenames1, read.table,header=T,fill=T)
id_wave=330
id_wave=as.character(id_wave)
track_data_wave=paste("/home/users/asawigr/reef/zad15/Dane/3_wysiew_liscie/",id_wave,"nm",sep="")
setwd(track_data_wave) 
filenames2 <- dir(pattern="arw",full.names=F)  
ilm=length(matrices) 
time=as.vector(matrices[[1]][2402:15603,1])
setwd(track_results) 
if (as.real(id_ch)==1){
write.table(time,"time.txt")}
L=13201
ptl=Sys.time()
print(ptl)
ilm
ilmall=ilm

#Tables for reference chromatograms prepereation - for all varieties 
REF1=NULL
REF2=NULL
IND_VARIETY=NULL
IND_LAB=NULL


##########The reference chromatogram selection inside varieties 
for (id_variety in 1:length(obk)){
#nrp-currently considering the variety for which we find reference chromatogram
#numkol - numbers from the data table set in order of loading files
#samnr - variety name according to the order in the data table
#nrtab - numbers successively loaded files on the selected variety
nrtab=NULL
Labels=NULL
countm=0
samnr=as.character(LAB[,2])
ilm=ilmall
for (ll1 in 1:ilm){
numkol=as.real(substr(as.character(matrices[[ll1]][1,2]),1,4))
numkol=which(LAB[,1]==numkol)
if (samnr[numkol]==obk[id_variety]){
nrtab=rbind(nrtab,ll1)
#Sample labels creation - reference chromatogram is not as a first, it is just on its place
Labels=rbind(Labels,LAB[numkol,])
#Samples from selected variety counting
countm=countm+1
}
}

print(countm)
ilm=countm


#matrices1 will include only selected variety for the wavelength 280nm
#matrices2 will include only selected variety for the wavelength 330nm
id_wave=280
id_wave=as.character(id_wave)
track_data_wave=paste("/home/users/asawigr/reef/zad15/Dane/3_wysiew_liscie/",id_wave,"nm",sep="")
setwd(track_data_wave) 
matrices1 <- lapply(filenames1[nrtab], read.table,header=T,fill=T)
id_wave=330
id_wave=as.character(id_wave)
track_data_wave=paste("/home/users/asawigr/reef/zad15/Dane/3_wysiew_liscie/",id_wave,"nm",sep="")
setwd(track_data_wave) 
matrices2 <- lapply(filenames2[nrtab], read.table,header=T,fill=T)

##########Data normalization
#All chromatograms are dividing by the mass of each chromatogram
masa=NULL
for (mi in 1:ilm){
masa=c(masa,matrices1[[mi]][1,3])
}
if (as.real(id_ch)==1){
setwd(track_results)
###Rys.a0.
pdf1=paste(as.character(id_variety),"_280_raw.pdf",sep="")
title1=paste("Raw data - ",obk[id_variety]," - 280 nm",sep="")
pdf(pdf1,width=10,height=6)
plot(time,as.vector(matrices1[[1]][2402:15603,2]),type="n",ylim=c(0,1),main=title1,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices1[[cc]][2402:15603,2]),col=cc)
}
dev.off()

###Rys.a0.
pdf1=paste(as.character(id_variety),"_330_raw.pdf",sep="")
title1=paste("Raw data - ",obk[id_variety]," - 330 nm",sep="")
pdf(pdf1,width=10,height=6)
plot(time,as.vector(matrices2[[1]][2402:15603,2]),type="n",ylim=c(0,1),main=title1,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices2[[cc]][2402:15603,2]),col=cc)
}
dev.off()
}


for (mi in 1:ilm){
matrices1[[mi]][2402:15603,2]=matrices1[[mi]][2402:15603,2]/masa[mi]
matrices2[[mi]][2402:15603,2]=matrices2[[mi]][2402:15603,2]/masa[mi]
}

if (as.real(id_ch)==1){
###Rys.a.
pdf1=paste(as.character(id_variety),"_280_after_normalization.pdf",sep="")
title1=paste(obk[id_variety]," after normalization - 280 nm",sep="")
pdf(pdf1,width=10,height=6)
plot(time,as.vector(matrices1[[1]][2402:15603,2]),type="n",ylim=c(0,3),main=title1,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices1[[cc]][2402:15603,2]),col=cc)
}
dev.off()

###Rys.a.
pdf1=paste(as.character(id_variety),"_330_after_normalization.pdf",sep="")
title1=paste(obk[id_variety]," after normalization - 330 nm",sep="")
pdf(pdf1,width=10,height=6)
plot(time,as.vector(matrices2[[1]][2402:15603,2]),type="n",ylim=c(0,3),main=title1,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices2[[cc]][2402:15603,2]),col=cc)
}
dev.off()
}

########Baseline removal based on differentiation:
for (mi in 1:ilm){
matrices1[[mi]][2402:15602,2]=diff(matrices1[[mi]][2402:15603,2])
matrices2[[mi]][2402:15602,2]=diff(matrices2[[mi]][2402:15603,2])
}

if (as.real(id_ch)==1){
###Rys.b.

pdf2a=paste(as.character(id_variety),"_280_before_cow_spr.pdf",sep="")
title2a=paste(obk[id_variety]," - before cow_spr - 280 nm",sep="")
pdf(pdf2a,width=10,height=6)
plot(time,as.vector(matrices1[[1]][2402:15603,2]),type="n",ylim=c(-0.06,0.06),main=title2a,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices1[[cc]][2402:15603,2]),col=cc)
}
dev.off()

###Rys.b.

pdf2a=paste(as.character(id_variety),"_330_before_cow_spr.pdf",sep="")
title2a=paste(obk[id_variety]," - before cow_spr - 330 nm",sep="")
pdf(pdf2a,width=10,height=6)
plot(time,as.vector(matrices2[[1]][2402:15603,2]),type="n",ylim=c(-0.06,0.06),main=title2a,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices2[[cc]][2402:15603,2]),col=cc)
}
dev.off()
}

########Similarity index (sim_ind) calculation for each biological sample for the two wavelengths separately
sim_ind=matrix(1,1,ilm)
for (aa in 1:ilm){
ref=matrices1[[aa]][2402:15602,2]
for (bb in 1:ilm){
sec=matrices1[[bb]][2402:15602,2]
sim_ind[1,aa]=sim_ind[1,aa]*abs(cor(ref, sec, method = c("pearson")))
}
}
sim_ind1=sim_ind

sim_ind=matrix(1,1,ilm)
for (aa in 1:ilm){
ref=matrices2[[aa]][2402:15602,2]
for (bb in 1:ilm){
sec=matrices2[[bb]][2402:15602,2]
sim_ind[1,aa]=sim_ind[1,aa]*abs(cor(ref, sec, method = c("pearson")))
}
}
sim_ind2=sim_ind

##########Similarity index (sim_ind) calculation for each biological sample as an average of the two wavelengths within the variety
for (i in 1:ilm){
sim_ind[i]=mean(c(sim_ind1[i],sim_ind2[i]))}
ind_ref=which.max(sim_ind)
print(ind_ref)
####
REF1=rbind(REF1,c(ind_ref,matrices1[[ind_ref]][2402:15602,2]))
REF2=rbind(REF2,c(ind_ref,matrices2[[ind_ref]][2402:15602,2]))
IND_VARIETY=rbind(IND_VARIETY,ind_ref)
IND_LAB=rbind(IND_LAB,Labels[ind_ref,])
}

setwd(track_results)
REFinfo=cbind(IND_VARIETY,IND_LAB)
write.table(REFinfo,"REF_info_rys_ident.txt")
if (as.real(id_ch)==1){
REFinfo=cbind(IND_VARIETY,IND_LAB)
write.table(REFinfo,"REF_info.txt")
write.table(REF1[,2:ncol(REF1)],"REF_280.txt")
write.table(REF2[,2:ncol(REF2)],"REF_330.txt")}

#write.csv(IND_VARIETY,"ref_ind_variety.csv")
#write.csv(IND_LAB,"ref_lab_variety.csv")

##########################################################################

########The reference chromatogram selection between varieties
ilm=nrow(REF1)
sim_ind=matrix(1,1,ilm)
for (aa in 1:ilm){
ref=REF1[aa,2:(L+1)]
for (bb in 1:ilm){
sec=REF1[bb,2:(L+1)]
sim_ind[1,aa]=sim_ind[1,aa]*abs(cor(ref, sec, method = c("pearson")))
}
}
sim_ind1=sim_ind

sim_ind=matrix(1,1,ilm)
for (aa in 1:ilm){
ref=REF2[aa,2:(L+1)]
for (bb in 1:ilm){
sec=REF2[bb,2:(L+1)]
sim_ind[1,aa]=sim_ind[1,aa]*abs(cor(ref, sec, method = c("pearson")))
}
}
sim_ind2=sim_ind


###########Similarity index (sim_ind) calculation for reference samples from each variety as the average of the two wavelength - one reference sample is defined amoung all samples
for (i in 1:ilm){
sim_ind[i]=mean(c(sim_ind1[i],sim_ind2[i]))}
ind_ref=which.max(sim_ind)
ind_ref=42
print("reference variety:")
print(ind_ref)
print(obk[ind_ref])
if (as.real(id_ch)==1){write.table(ind_ref,"id_line_ref.txt")}  #saving in the file the number of the reference line / variety from the list of lines and varieties obk


#ref1, ref2 - reference chromatograms for all chromatograms ref1 - for the wavelength 280nm, ref1 - for the wavelength 330nm 
ref1=REF1[ind_ref,2:(L+1)]
ref2=REF2[ind_ref,2:(L+1)]
###########


time=time[1:L]

setwd('/home/users/asawigr/reef/zad15/3_wysiew_liscie/warping') 


pdf("42_ref_org_po_cow_MCam068_280.pdf",width=10,height=6)
plot(time,REF1[ind_ref,2:(L+1)]),type="l",main="MCam068_280", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()

pdf("42_ref_org_po_cow_MCam068_280_fragment1.pdf",width=10,height=6)
plot(time,REF1[ind_ref,2:(L+1)]),type="l",xlim=c(2.5,4.5),main="MCam068_280", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()

pdf("42_ref_org_po_cow_MCam068_280_fragment2.pdf",width=10,height=6)
plot(time,REF1[ind_ref,2:(L+1)]),type="l",xlim=c(4.5,6.5),main="MCam068_280", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()

pdf("42_ref_org_po_cow_MCam068_280_fragment3.pdf",width=10,height=6)
plot(time,REF1[ind_ref,2:(L+1)]),type="l",xlim=c(6.5,8.5),main="MCam068_280", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()

pdf("42_ref_org_po_cow_MCam068_280_fragment4.pdf",width=10,height=6)
plot(time,REF1[ind_ref,2:(L+1)]),type="l",xlim=c(8.5,10.5),main="MCam068_280", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()

pdf("42_ref_org_po_cow_MCam068_280_fragment5.pdf",width=10,height=6)
plot(time,REF1[ind_ref,2:(L+1)]),type="l",xlim=c(10.5,13),main="MCam068_280", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()



pdf("42_ref_org_po_cow_MCam068_330.pdf",width=10,height=6)
plot(time,REF2[ind_ref,2:(L+1)]),type="l",main="MCam068_330", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()

pdf("42_ref_org_po_cow_MCam068_330_fragment1.pdf",width=10,height=6)
plot(time,REF2[ind_ref,2:(L+1)]),type="l",xlim=c(2.5,4.5),main="MCam068_330", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()

pdf("42_ref_org_po_cow_MCam068_330_fragment2.pdf",width=10,height=6)
plot(time,REF2[ind_ref,2:(L+1)]),type="l",xlim=c(4.5,6.5),main="MCam068_330", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()

pdf("42_ref_org_po_cow_MCam068_330_fragment3.pdf",width=10,height=6)
plot(time,REF2[ind_ref,2:(L+1)]),type="l",xlim=c(6.5,8.5),main="MCam068_330", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()

pdf("42_ref_org_po_cow_MCam068_330_fragment4.pdf",width=10,height=6)
plot(time,REF2[ind_ref,2:(L+1)]),type="l",xlim=c(8.5,10.5),main="MCam068_330", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()

pdf("42_ref_org_po_cow_MCam068_330_fragment5.pdf",width=10,height=6)
plot(time,REF2[ind_ref,2:(L+1)]),type="l",xlim=c(10.5,13),main="MCam068_330", xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]",col=1,axes=FALSE)
axis(1,seq(2,11,0.1))
dev.off()