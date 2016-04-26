#The reference chromatogram selection and reference chromatograms for each variety after COW at the grid export
#1. The reference chromatogram selection inside varieties
#1.1. Data normalization
#1.2. Baseline removal based on differentiation
#1.3. Similarity index (sim_ind) calculation for each biological sample for the two wavelengths separately
#1.4. Similarity index (sim_ind) calculation for each biological sample as an average of the two wavelengths within the variety
#2. The reference chromatogram selection between varieties
#2.1. Similarity index (sim_ind) calculation for reference samples from each variety as the average of the two wavelength - one reference sample is defined amoung all samples
#3. Reference chromatograms for each variety after COW at the grid export
#3.1. Warping for reference chromatograms for varieties separately for both wavelength 
#3.2. Reference chromatograms for varieties after warping for normalized and diffarenced data saving in the files 001_cow_grid_wave280.txt, 002_cow_grid_wave280.txt,...

mv=c(132,176,220,275)
tv=c(10,12,15,18)

track_data_general=paste("/home/users/asawigr/reef/zad15/Dane/3_wysiew_liscie/",sep="")
track_results=paste("/home/users/asawigr/reef/zad15/3_wysiew_liscie/warping",sep="")
setwd(track_data_general)

#Table of observation import
LAB=read.csv("tabela_obserwacji_zad15_doœw3.csv")
obk=unique(LAB[,2])

#Data import
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


##########1. The reference chromatogram selection inside varieties 
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


##########1.1. Data normalization
#All chromatograms are dividing by the mass of each chromatogram
masa=NULL
for (mi in 1:ilm){
masa=c(masa,matrices1[[mi]][1,3])
}
if (as.real(id_ch)==1){
setwd(track_results)
#Rys.a0.
pdf1=paste(as.character(id_variety),"_280_raw.pdf",sep="")
title1=paste("Raw data - ",obk[id_variety]," - 280 nm",sep="")
pdf(pdf1,width=10,height=6)
plot(time,as.vector(matrices1[[1]][2402:15603,2]),type="n",ylim=c(0,1),main=title1,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices1[[cc]][2402:15603,2]),col=cc)
}
dev.off()

#Rys.a0.
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
#Rys.a.
pdf1=paste(as.character(id_variety),"_280_after_normalization.pdf",sep="")
title1=paste(obk[id_variety]," after normalization - 280 nm",sep="")
pdf(pdf1,width=10,height=6)
plot(time,as.vector(matrices1[[1]][2402:15603,2]),type="n",ylim=c(0,3),main=title1,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices1[[cc]][2402:15603,2]),col=cc)
}
dev.off()

#Rys.a.
pdf1=paste(as.character(id_variety),"_330_after_normalization.pdf",sep="")
title1=paste(obk[id_variety]," after normalization - 330 nm",sep="")
pdf(pdf1,width=10,height=6)
plot(time,as.vector(matrices2[[1]][2402:15603,2]),type="n",ylim=c(0,3),main=title1,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices2[[cc]][2402:15603,2]),col=cc)
}
dev.off()
}


##########1.2. Baseline removal based on differentiation
for (mi in 1:ilm){
matrices1[[mi]][2402:15602,2]=diff(matrices1[[mi]][2402:15603,2])
matrices2[[mi]][2402:15602,2]=diff(matrices2[[mi]][2402:15603,2])
}

if (as.real(id_ch)==1){
#Rys.b.
pdf2a=paste(as.character(id_variety),"_280_before_cow_spr.pdf",sep="")
title2a=paste(obk[id_variety]," - before cow_spr - 280 nm",sep="")
pdf(pdf2a,width=10,height=6)
plot(time,as.vector(matrices1[[1]][2402:15603,2]),type="n",ylim=c(-0.06,0.06),main=title2a,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices1[[cc]][2402:15603,2]),col=cc)
}
dev.off()

#Rys.b.
pdf2a=paste(as.character(id_variety),"_330_before_cow_spr.pdf",sep="")
title2a=paste(obk[id_variety]," - before cow_spr - 330 nm",sep="")
pdf(pdf2a,width=10,height=6)
plot(time,as.vector(matrices2[[1]][2402:15603,2]),type="n",ylim=c(-0.06,0.06),main=title2a,xlab="Time [min]",ylab="Intensity of absorbance [A.U.]",col=2)
for (cc in 1:ilm){
lines(time,as.vector(matrices2[[cc]][2402:15603,2]),col=cc)
}
dev.off()
}


##########1.3. Similarity index (sim_ind) calculation for each biological sample for the two wavelengths separately
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


##########1.4. Similarity index (sim_ind) calculation for each biological sample as an average of the two wavelengths within the variety
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
if (as.real(id_ch)==1){
REFinfo=cbind(IND_VARIETY,IND_LAB)
write.table(REFinfo,"REF_info.txt")
write.table(REF1[,2:ncol(REF1)],"REF_280.txt")
write.table(REF2[,2:ncol(REF2)],"REF_330.txt")}

#write.csv(IND_VARIETY,"ref_ind_variety.csv")
#write.csv(IND_LAB,"ref_lab_variety.csv")

##########################################################################

##########2. The reference chromatogram selection between varieties
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


##########2.1. Similarity index (sim_ind) calculation for reference samples from each variety as the average of the two wavelength - one reference sample is defined amoung all samples
for (i in 1:ilm){
sim_ind[i]=mean(c(sim_ind1[i],sim_ind2[i]))}
ind_ref=which.max(sim_ind)
print("reference variety:")
print(ind_ref)
print(obk[ind_ref])
if (as.real(id_ch)==1){write.table(ind_ref,"id_line_ref.txt")} #saving in the file the number of the reference line / variety from the list of lines and varieties obk


#ref1, ref2 - reference chromatograms for all chromatograms ref1 - for the wavelength 280nm, ref1 - for the wavelength 330nm 
ref1=REF1[ind_ref,2:(L+1)]
ref2=REF2[ind_ref,2:(L+1)]

if (as.real(id_ch)==1){
#Plot before warping
setwd(track_results)
pdf("ref_all_280.pdf",width=10,height=6)
plot(time[1:L],ref,ylim=c(-0.06,0.06),type="n",main="Reference chromatograms - 280nm",col=1)
for (cc in 1:ilm){
lines(time[1:L],REF1[cc,2:(L+1)],col=cc)
}
dev.off()

pdf("ref_all_330.pdf",width=10,height=6)
plot(time[1:L],ref,ylim=c(-0.06,0.06),type="n",main="Reference chromatograms - 330nm",col=1)
for (cc in 1:ilm){
#if (cc!=ind_ref){
lines(time[1:L],REF2[cc,2:(L+1)],col=cc)
}
#}
dev.off()
}


time=time[1:L]

###################################################

##########3. Reference chromatograms for each variety after COW at the grid export

#Function cow_opt_grid for the pair of chromatograms

fun_cow_opt_grid=function(L,time,ref,sec,lsec,mv,tv,id_wave){
d1=ref  
d2=sec

#M=cbind(d1, d2)
#matplot(time,M, type="l", main = "Baseline removal based on differentiation", col=1:2, lty=1)
#legend("topright",paste(c("ref","sec")), col=1:2, lty=1)

#################
#Optimization

Lp=(L-1)
R=NULL
#ptl=Sys.time()
#print(ptl)

for (countm in 1:length(mv)){
for (countt in 1:length(tv)){ 
m=mv[countm]
t=tv[countt]
print(m)
print(t)
ptl=Sys.time()
print(ptl)
if (m>(t+3)) {
N=Lp/m

###################
#COW
T=d1
P=d2

###############################################
time_o=c(1:L)
delta=0
LT=m*N
fsum=0
F=matrix(-Inf,(N+1),(LT+1))
U=F
F[N+1,N*m+1]=0
	
for (i in N:1) {
s1=max((i-1)*(m+delta-t),LT-(N-i+1)*(m+delta+t))+1
e1=min((i-1)*(m+delta+t),LT-(N-i+1)*(m+delta-t))+1
for (x in s1:e1) { for (u in (delta-t):(delta+t)){
#for different u different F(i,x) is counting and max is chosen

#linear interpolation
#y=y0+(x-x0)(y1-y0)/(x1-x0)
#s2, e2 - the begining and the end of interval in P'
#s, e - the begining and the end of interval in P
s2=x

k=x+m+u
s=m*i-m+1
e=m*(i+1)-m+1
if (k>=(m*N+1))  ###break the loop with u
k=m*N+1

e2=k

p=rep(0,e2-s2+1)

for (j in 0:(e2-s2)) {
p[j+1]=j*(e-s)/(e2-s2)+s
}

#cat("New time positions")
#print(p) #New time positions

P1=approx(time_o,P, xout=p, method="linear")$y
#cat("Warped P, czyli P'")
#print(P1)

##############

#Pearson correlation coefficient (c) counting in the interval (x;x+m+u)
#c=cor(P1[x:(x+m+u)], T[x:(x+m+u)], method = c("pearson"))
len=length(P1)
if ((P1==rep(0,len)) && (T[x:k]==rep(0,len)))
c=1 else
c=cor(P1, T[x:k], method = c("pearson"))


#cat("i:")
#print(i)
#cat("x:")
#print(x)
#cat("k:")
#print(k)
#cat("P1:")
#print(P1)
#cat("T:")
#print(T[x:k])

#print(c)

fsum=F[(i+1),k]+c
#if (i==N) fsum=c
#fsum=fsum+c
#print(fsum)

if (fsum>F[i,x]) {
F[i,x]=fsum
U[i,x]=u
}
}
}              
}
#print(F)
#print(U)

x[1]=1
for (i in 1:N) {
u[i]=U[i,x[i]]
x[i+1]=x[i]+m+u[i]
}
#x=optimal x*
#print(x)

P1=rep(0,N*m+1)
##P' recovery
poprz=1
for (i in 1:N){
s2=x[i]

s=m*i-m+1
e=m*(i+1)-m+1
e2=x[i+1]

p=rep(0,e2-s2+1)
for (j in 0:(e2-s2)) {
p[j+1]=j*(e-s)/(e2-s2)+s
}
#cat("New time positions")
#print(p) #New time positions
P2=approx(time_o,P, xout=p, method="linear")$y
#cat("P2: ")
#print(P2)
#print(poprz)
P1[(poprz):(poprz+e2-s2)]=P2
poprz=poprz+e2-s2
#cat("P1: ")
#print(P1)
}

row=cbind(m,t,t(P1))
R=rbind(R,row)

#Plot - results
#if (t==1){
#M=cbind(T, P, P1)
#matplot(time,M, type="l", main = "Aligment using COW",xlab = "Czas [min]", ylab= "Intensywnosc absorbcji [A.U.]", col=1:3, lty=1)
#legend("topright",paste(c("ref (T)","sec przed COW (P)","sec po zast. COW (P')")), col=1:3, lty=1)
#}
##
}
}
}
text0=as.character(lsec)
if (nchar(text0)==3){
text1=paste("gridref_",as.character(id_wave),"_",text0,sep="")}
if (nchar(text0)==2){
text1=paste("gridref_",as.character(id_wave),"_0",text0,sep="")}
if (nchar(text0)==1){
text1=paste("gridref_",as.character(id_wave),"_00",text0,sep="")}
text2=paste(text1,".txt",sep="")
write.table(R,text2)


print("The correct ending")
ptl=Sys.time()
print(ptl)
} #the end of funtion fun_cow_opt_grid


##########3.1. Warping for reference chromatograms for varieties separately for both wavelength 

###cow_opt_grid for reference chromatograms
id_wave=as.real(id_wave_ref) #loading condidered wavelength
if (id_wave==280){ref=ref1
REF=REF1} else {ref=ref2
REF=REF2}

for (bb in as.real(id_ch):as.real(id_ch)){
if (bb!=ind_ref){
sec=REF[bb,2:(L+1)]
lsec=bb
#function cow_opt for the pair of chromatograms

fun_cow_opt_grid(L,time,ref,sec,lsec,mv,tv,id_wave)
}
}

##########3.2. Reference chromatograms for varieties after warping for normalized and diffarenced data saving in the files 001_cow_grid_wave280.txt, 002_cow_grid_wave280.txt,...

text0=as.character(ind_ref)
if (nchar(text0)==3){
text1=paste("gridref_",as.character(id_wave),"_",text0,sep="")}
if (nchar(text0)==2){
text1=paste("gridref_",as.character(id_wave),"_0",text0,sep="")}
if (nchar(text0)==1){
text1=paste("gridref_",as.character(id_wave),"_00",text0,sep="")}
text2=paste(text1,".txt",sep="")
write.table("0",text2)