#Common peaks for all chromatograms detection and quantitation

#1. Second derivative and its smoothing
#2. Peak detection on the basis of the sign of the differented chromatograms
#3. The chromatograms integration over peaks for all peaks and all samples


L=13201

track_results=paste("/home/users/asawigr/reef/zad15/3_wysiew_liscie/warping",sep="")

setwd(track_results)

variety_ref=read.table("id_line_ref.txt") #the reference variety number loading
ind_ref=as.real(variety_ref)
ind_ref

time=as.matrix(read.table("time.txt"))
time=time[1:L] 

ref_info=read.table("REF_info.txt")
ref_280=ref_info[,1]
ref_330=ref_280
#new X after warping
X=NULL
Z=NULL
filenames <- dir(pattern="PCAmatrix_grid_280",full.names=F)  
U <- lapply(filenames, read.table,header=T,fill=T)

for (i in 1:length(filenames)){
w1=U[[i]][1,]
if (ref_280[i]!=1){
U[[i]][1:(ref_280[i]-1),]=U[[i]][2:ref_280[i],]
U[[i]][ref_280[i],]=w1}
X=rbind(X,as.matrix(U[[i]]))
}

filenames <- dir(pattern="PCAmatrix_grid_330",full.names=F)  
Y <- lapply(filenames, read.table,header=T,fill=T)

for (i in 1:length(filenames)){
w1=Y[[i]][1,]
if (ref_280[i]!=1){
Y[[i]][1:(ref_280[i]-1),]=Y[[i]][2:ref_280[i],]
Y[[i]][ref_280[i],]=w1}
Z=rbind(Z,as.matrix(Y[[i]]))
}


#######Csv file preparation
#labels and ind_ref import
#conversion to pre-order in the labels, the ref will not be the first, to order samples was the same for both wavelengths
#analogous conversion of the results, that is, in Xn ref goes to his place
wave=280
LABELS2=NULL

filenames <- dir(pattern=paste("labels_samp_",as.character(wave),sep=""),full.names=F) 
print(filenames) 
LABELS <- lapply(filenames, read.table,header=T,fill=T)

for (i in 1:length(filenames)){
LABELS2=rbind(LABELS2,as.matrix(LABELS[[i]][,2:5]))
}

#######Other terms for drought treatments, time point of observations and replications:

for (i in 1:nrow(LABELS2)){
if (as.character(LABELS2[i,4])=="S1"){LABELS2[i,4]="I"}
if (as.character(LABELS2[i,4])=="S2"){LABELS2[i,4]="II"}
}


###########Function - the choice of the peaks and chromatograms integration over peaks for X (wavelength 280nm) and Y (wavelength 330nm)
fun_peaks=function(X,wave,time,L){
wave=as.character(wave)
##########1.a. Second derivative
#Baseline removal based on differentiation
Diff2b=NULL 
for (a in 1:nrow(X)){
Xb=diff(X[a,1:L]) 
Diff2b=rbind(Diff2b,Xb)
}
#pdf(paste("before_smoothing_",wave,".pdf",sep=""),width=10,height=6)
#part=t(Diff2b)
#matplot(time[1:(L-1)],part[,1:10],type="l",ylim=c(-0.002,0.002),xlim=c(4.5,6),main="Second derivative (after cow) without smoothing")
#dev.off()

##########1.b. Second derivative smoothing with n=520
Diff2=NULL
for (a in 1:nrow(X)){
sp=smooth.spline(time[1:(L-1)], Diff2b[a,1:(L-1)],n=520)
Diff2=rbind(Diff2,sp$y)
}
#pdf(paste("after_smoothing_",wave,".pdf",sep=""),width=10,height=6)
#part=t(Diff2)
#matplot(time[1:(L-1)],part[,1:10],type="l",ylim=c(-0.002,0.002),xlim=c(4.5,6),main="Second derivative (after cow) with smoothing")
#dev.off()


###########2. Peak detection on the basis of the sign of the differented chromatograms

#########
#Peak selection (where second derivative negative)
Signs=sign(Diff2)
c=rep(1:1,times=nrow(X))
SignsPlus=cbind(c,Signs,c)
Peaksall=NULL
elPeaks=matrix(0,1,nrow(X))
for (j in 1:nrow(X)){
print(j)
PeaksL=NULL
PeaksR=NULL
Peaks=NULL
#Peaks - left and right bound of interval and time
for (i in 2:(L+1)){
if ((SignsPlus[j,i]==-1)&&(SignsPlus[j,(i-1)]==1)){
PeaksL=rbind(PeaksL,c(i,time[i]))
}
if ((SignsPlus[j,i]==1)&&(SignsPlus[j,(i-1)]==-1)){
PeaksR=rbind(PeaksR,c((i-1),time[(i-1)]))
}
}
Peaks=cbind(PeaksL[,1],PeaksR[,1],PeaksL[,2],PeaksR[,2])

#################The peaks of length 1 deletion
peak_length=1
Dpeaks=NULL
for (yy in 1:nrow(Peaks)){
if (as.real(Peaks[yy,2])-as.real(Peaks[yy,1])>peak_length){
Dpeaks=rbind(Dpeaks,Peaks[yy,])}
}
Peaks=Dpeaks
######################

##########Tolerances to remove the peaks
#if max element in the peak < tol, than this peak is eliminated from the list of peaks
tol=0.0005
####Irrelevant peaks removing
Peaks2=NULL
if (Peaks[nrow(Peaks),2]==L){
Peaks[nrow(Peaks),2]=(L-1)
}
for (i in 1:nrow(Peaks)){
mind=max(X[j,(Peaks[i,1]:Peaks[i,2])])
if (abs(mind)>tol){
Peaks2=rbind(Peaks2,Peaks[i,])
}
}

Peaks=Peaks2
Peaksall=rbind(Peaksall,Peaks)
if (length(Peaks)!=0){
elPeaks[j]=nrow(Peaks2)}
} #end for j

#####################################################
# All - the vector of zero elements, for each element 1 is added if there is a peak in a single chromatogram;
# the result is a vector with element, which are the numbers of chromagrams (for all chromatograms, for all varieties), for which at a given time point there is a peak
# i.e. those that have a negative II-nd difference (II-nd derivative)
All=matrix(0,1,(L-1))
for (i in 1:nrow(Peaksall)){
for (j in Peaksall[i,1]:Peaksall[i,2]){ 
All[j]=All[j]+1
}
}
########

#################
#Peaks selection
#Intervals from All are selected, which all element are greater than mch = 10, and for too long intervals (max_peak = 0.1) mch = 50 is used
#which means that: those peaks are selected that occur for at least mch = 10 chromatograms
#and too long peaks if it can be divided with mch = 50 they are divided, if not they stay not changed.
AllPlus=cbind(1,All,1)
PeaksLA=NULL
PeaksRA=NULL
PeaksA=NULL
mch=10 #the minimum number of chromatograms in which a peak occure to be recognized
#Peaks - left and right bound of interval and time
for (i in 2:(L+1)){
if ((AllPlus[i]>=mch)&&(AllPlus[(i-1)]<mch)){
PeaksLA=rbind(PeaksLA,c(i,time[i]))
}
if ((AllPlus[i]<mch)&&(AllPlus[(i-1)]>=mch)){
PeaksRA=rbind(PeaksRA,c((i-1),time[(i-1)]))
}
}
Peakscommon=cbind(PeaksLA[,1],PeaksRA[,1],PeaksLA[,2],PeaksRA[,2])
Peakscommon

# tol - minimum peak height (or more precisely the difference after warping) in a single chromatogram
# mch - the minimum number of chromatograms that make up the common peak
# max_peak - # maximum width of the common peak
# if max_peak is exceeded, it will attempt to divide the peak
# for more peaks, i.e. in such intervals mch increasing to 30

Peakscommonold=Peakscommon
max_peak=0.1  
mch=50
sump=0
for (i in 1:nrow(Peakscommonold)){
PeaksLAp=NULL
PeaksRAp=NULL
if ((as.real(Peakscommonold[i,4])-as.real(Peakscommonold[i,3]))>max_peak){
for (ii in (as.real(Peakscommonold[i,1])+1):(as.real(Peakscommonold[i,2])+1)){
if ((AllPlus[ii]>=mch)&&(AllPlus[ii-1]<mch)){
PeaksLAp=rbind(PeaksLAp,c(ii,time[ii]))
}
if ((AllPlus[ii]<mch)&&(AllPlus[ii-1]>=mch)){
PeaksRAp=rbind(PeaksRAp,c((ii-1),time[(ii-1)]))
}
} #for ii
if ((length(PeaksRAp)!=0)||(length(PeaksLAp)!=0)){
if (nrow(PeaksRAp)==nrow(PeaksLAp)){
Peakscommonpoprz=Peakscommon
indp1=i+sump
indp2=i+sump-1+nrow(PeaksRAp)
Peakscommon[indp1:indp2,1]=PeaksLAp[,1]
Peakscommon[indp1:indp2,2]=PeaksRAp[,1]
Peakscommon[indp1:indp2,3]=PeaksLAp[,2]
Peakscommon[indp1:indp2,4]=PeaksRAp[,2]
Peakscommon=rbind(Peakscommon,matrix(0,(nrow(PeaksLAp)-1),4))
if (i<nrow(Peakscommonold)) {Peakscommon[(indp2+1):(nrow(Peakscommonpoprz)-1+nrow(PeaksLAp)),]=Peakscommonold[(i+1):nrow(Peakscommonold),]}
sump=sump+nrow(PeaksLAp)-1
} #if nrow()<>NULL
} #if PeaksRAp<>NULL and PeaksLAp<>NULL
} #if peak longer than max_peak
} #for i 


##################The peaks of length <=10 removing
peak_length=10
Dlpeaks=NULL
for (yy in 1:nrow(Peakscommon)){
if (as.real(Peakscommon[yy,2])-as.real(Peakscommon[yy,1])>peak_length){
Dlpeaks=rbind(Dlpeaks,Peakscommon[yy,])}
}
Peakscommon=Dlpeaks
######################


Peakscommon1<<-Peakscommon



##########3. The chromatograms integration over peaks for all peaks and all samples
# Labels - labels for samples in X
# The labels for the peaks are in the file lista_pików...txt
# TC - table with samples and peaks filled with integrals
TC=matrix(0,nrow(X),nrow(Peakscommon))
Y=t(X)
Y=abs(Y)  
x1=as.real(time[1:L])
for (i in 1:nrow(X)){
y1=as.real(Y[,i])
for (j in 1:nrow(Peakscommon)){
x=x1[Peakscommon[j,1]:Peakscommon[j,2]]
y=y1[Peakscommon[j,1]:Peakscommon[j,2]]
suma=sum(y)
TC[i,j]=suma
}
}

#################
text0=paste("peak_detection_",wave,"_10_tol0_0005_max_peak0_1_50.pdf",sep="")
pdf(text0,width=10,height=6)
matplot(time[1:L],t(X),ylim=c(-0.06,0.06),type="l",main=paste("Dane po cow z pikami - ",wave,sep=""))
abline(v=Peakscommon[,3])
abline(v=Peakscommon[,4],col="red")
dev.off()

text1=paste("peak_detection_",wave,"part1_10_tol0_0005_max_peak0_1_50.pdf",sep="")
pdf(text1,width=10,height=6)
matplot(time[1:L],t(X),ylim=c(-0.06,0.06),xlim=c(2.5,4.5),type="l",main=paste("Dane po cow z pikami",wave," - part",sep=""))
abline(v=Peakscommon[,3])
abline(v=Peakscommon[,4],col="red")
dev.off()

text2=paste("peak_detection_",wave,"part2_10_tol0_0005_max_peak0_1_50.pdf",sep="")
pdf(text2,width=10,height=6)
matplot(time[1:L],t(X),ylim=c(-0.06,0.06),xlim=c(4.5,6.5),type="l",main=paste("Dane po cow z pikami",wave," - part",sep=""))
abline(v=Peakscommon[,3])
abline(v=Peakscommon[,4],col="red")
dev.off()

text3=paste("peak_detection_",wave,"part3_10_tol0_0005_max_peak0_1_50.pdf",sep="")
pdf(text3,width=10,height=6)
matplot(time[1:L],t(X),ylim=c(-0.06,0.06),xlim=c(6.5,8.5),type="l",main=paste("Dane po cow z pikami",wave," - part",sep=""))
abline(v=Peakscommon[,3])
abline(v=Peakscommon[,4],col="red")
dev.off()

text4=paste("peak_detection_",wave,"part4_10_tol0_0005_max_peak0_1_50.pdf",sep="")
pdf(text4,width=10,height=6)
matplot(time[1:L],t(X),ylim=c(-0.06,0.06),xlim=c(8.5,10.5),type="l",main=paste("Dane po cow z pikami",wave," - part",sep=""))
abline(v=Peakscommon[,3])
abline(v=Peakscommon[,4],col="red")
dev.off()

text5=paste("peak_detection_",wave,"part5_10_tol0_0005_max_peak0_1_50.pdf",sep="")
pdf(text5,width=10,height=6)
matplot(time[1:L],t(X),ylim=c(-0.06,0.06),xlim=c(10.5,13),type="l",main=paste("Dane po cow z pikami",wave," - part",sep=""))
abline(v=Peakscommon[,3])
abline(v=Peakscommon[,4],col="red")
dev.off()


TC
}#end of function fun_peaks

#funtion fun_peaks for wavelength 280nm i 330nm
TC1=fun_peaks(X,280,time,L)
peaks1=Peakscommon1
write.csv(TC1,"TC1spr.csv",row.names=FALSE)
Peakscommon1=cbind(matrix(280,nrow(peaks1),1),peaks1)
Peaks_list=Peakscommon1
Peakscommom1=NULL

TC2=fun_peaks(Z,330,time,L)
peaks2=Peakscommon1
Peakscommon1=cbind(matrix(330,nrow(peaks2),1),peaks2)
Peaks_list=rbind(Peaks_list,Peakscommon1)

Peaks_list=cbind(c(1:nrow(Peaks_list)),Peaks_list)
l=c("numpik","wave","numtmin","numtmax","tmin","tmax")
Peaks_list=rbind(l,Peaks_list)
write.csv(Peaks_list,"peaks.csv",row.names=FALSE)

TC=cbind(c(1:nrow(TC1)),LABELS2,TC1,TC2)

#TC=cbind(c(1:nrow(X)),TC)
p=NULL
for (i in 1:(nrow(Peaks_list)-1)){
p=cbind(p,paste("v",i,sep=""))}
l=c("numproby","obiekt!","termin!","lisc!","wariant!",p)
TC=rbind(as.matrix(t(l)),as.matrix(TC))
write.csv(TC,"dane.csv",row.names=FALSE)


l=NULL
for (i in 1:(nrow(Peaks_list))){
l=rbind(l,paste("TO:0000281_UPLC_",Peaks_list[i,2],"_",Peaks_list[i,5],"_",Peaks_list[i,6],sep=""))
}
Peaks_list=cbind(Peaks_list, l)

u=format(as.real(Peaks_list[2:nrow(Peaks_list),5]),digits=7)
v=format(as.real(Peaks_list[2:nrow(Peaks_list),6]),digits=7)
ind1=which(as.real(Peaks_list[2:nrow(Peaks_list),5])>9.999999)
ind2=which(as.real(Peaks_list[2:nrow(Peaks_list),6])>9.999999)
ind1=ind1
ind2=ind2

u[ind1]=format(as.real(Peaks_list[ind1+rep(1,length(ind1)),5]),digits=7)
v[ind2]=format(as.real(Peaks_list[ind2+rep(1,length(ind2)),6]),digits=7)

l=NULL
for (i in 2:(nrow(Peaks_list))){
ptext=paste("TO:0000281_UPLC_",Peaks_list[i,2],"_",u[i-1],"_",v[i-1],sep="")
ptext1=substr(ptext,1,20)
if (any((i-1)==ind1)) {ptext2=substr(ptext,21,30)} else {ptext2=substr(ptext,22,30)}
if (any((i-1)==ind2)) {ptext3=substr(ptext,31,39)} else {ptext3=substr(ptext,32,39)}
ptext=paste(ptext1,ptext2,ptext3,sep="")
l=rbind(l,ptext)
}
l=rbind("TO:0000281_UPLC_wave_tmin_tmax",l)
Peaks_list=cbind(Peaks_list, l)


write.csv(Peaks_list,"peaks_plus.csv",row.names=FALSE)

######To Database Germinate
l2=c("PhenotypeNo","Name","ShortName","Units")
Pht=cbind(c(1:(nrow(Peaks_list)-1)),Peaks_list[2:nrow(Peaks_list),8],rep("*",(nrow(Peaks_list)-1)),rep("bw",(nrow(Peaks_list)-1)))
Pht=rbind(l2,Pht)
write.csv(Pht,"Zad15_dosw1_phenotypes.csv",row.names=FALSE)

l3=c("Lp","Experiment","Task","Organ","Drought_treatment","Time_point","Germinate","Replicate",l[2:length(l)])
rc=NULL
for (i in 2:(nrow(TC))){
rc=cbind(rc,paste("R0",TC[i,4],sep=""))
}
Phd=cbind(c(1:(nrow(TC)-1)),rep("IGR_2011_1",(nrow(TC)-1)),rep("Z15",(nrow(TC)-1)),rep("Leaf",(nrow(TC)-1)),TC[2:nrow(TC),5],TC[2:nrow(TC),3],TC[2:nrow(TC),2],t(rc),TC[2:nrow(TC),6:ncol(TC)])
Phd=rbind(l3,Phd)
write.csv(Phd,"Zad15_dosw1_phenotypedata.csv",row.names=FALSE)


#pdf("test_5_sample_280.pdf",width=10,height=6)
#matplot(time[1:L],t(X[1:5,]),ylim=c(-0.02,0.02),xlim=c(2,4),type="l",main="Data after cow with peaks")
#abline(v=peaks1[,3])
#abline(v=peaks1[,4],col="red")
#dev.off()
