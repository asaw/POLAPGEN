#Peak division
#requires pvclust package

library(pvclust)

setwd("/home/users/asawigr/reef/zad15/3_wysiew_liscie/warping")
wave=330
Xall=read.csv(paste("X_",wave,".csv",sep=""))

peaks=read.csv("peaks.csv")
peaks=as.matrix(peaks[2:nrow(peaks),])
peaks=peaks[which(as.double(peaks[,2])==wave),]
peakst=peaks

table=read.csv("Zad15_dosw1_phenotypedata.csv")

if (wave==280){table=table[,1:(8+nrow(peaks))]}
if (wave==330){ncolt=ncol(table)
table=cbind(table[,1:8],table[,(ncolt-nrow(peaks)+1):ncolt])}


time=as.matrix(read.table("time.txt"))
L=ncol(Xall)


#The length of the signal calculation and too long of its indication, that is, longer than 0.1 min in the last column of peaks
peaks=cbind(peaks,(as.double(peaks[,6])-as.double(peaks[,5])))


fun_peaks=function(X,wave,time,L,longs){
wave=as.character(wave)
##Drugie ró¿nicowanie X
#Baseline removal based on differentiation
Diff2b=NULL 
for (a in 1:nrow(X)){
Xb=diff(as.double(X[a,1:L]))
Diff2b=rbind(Diff2b,Xb)
}


##########
#Second derivative smoothing with n=520
Diff2=NULL
for (a in 1:nrow(X)){
sp=smooth.spline(time[1:(L-1)], Diff2b[a,1:(L-1)],n=520)
Diff2=rbind(Diff2,sp$y)
}

############Peak detection on the basis of the sign of the differented chromatograms

#########
#Peak selection (where second derivative negative)
Signs=sign(Diff2)
c=rep(1:1,times=nrow(X))
SignsPlus=cbind(c,Signs,c)
Peaksall=NULL
elPeaks=matrix(0,1,nrow(X))
for (j in 1:nrow(X)){
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
if (as.double(Peaks[yy,2])-as.double(Peaks[yy,1])>peak_length){
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


if (!is.null(Peakscommon)){
# tol - minimum peak height (or more precisely the difference after warping) in a single chromatogram
# mch - the minimum number of chromatograms that make up the common peak
# max_peak - # maximum width of the common peak
# if max_peak is exceeded, it will attempt to divide the peak
# for more peaks, i.e. in such intervals mch increasing to 30
Peakscommonold=Peakscommon
max_peak=0.1  
mch=50
#If group consist of <=50 elements, than mch=10
if (nrow(X)<=50) {mch=10}
sump=0
for (i in 1:nrow(Peakscommonold)){
#print(i)
PeaksLAp=NULL
PeaksRAp=NULL
if ((as.double(Peakscommonold[i,4])-as.double(Peakscommonold[i,3]))>max_peak){
for (ii in (as.double(Peakscommonold[i,1])+1):(as.double(Peakscommonold[i,2])+1)){
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
if (indp2>nrow(Peakscommon)) {Peakscommon=rbind(Peakscommon,matrix(0,indp2-nrow(Peakscommon),4))}
Peakscommon[indp1:indp2,1]=PeaksLAp[,1]
Peakscommon[indp1:indp2,2]=PeaksRAp[,1]
Peakscommon[indp1:indp2,3]=PeaksLAp[,2]
Peakscommon[indp1:indp2,4]=PeaksRAp[,2]
Peakscommon=rbind(Peakscommon,matrix(0,(nrow(PeaksLAp)-1),4))

if (i<nrow(Peakscommonold)) {if ((nrow(Peakscommonpoprz)-1+nrow(PeaksLAp))>nrow(Peakscommon)) {Peakscommon=rbind(Peakscommon,matrix(0,nrow(Peakscommonpoprz)-1+nrow(PeaksLAp)-nrow(Peakscommon),4))}
Peakscommon[(indp2+1):(nrow(Peakscommonpoprz)-1+nrow(PeaksLAp)),]=Peakscommonold[(i+1):nrow(Peakscommonold),]}
sump=sump+nrow(PeaksLAp)-1
} #if nrow()<>NULL
} #if PeaksRAp<>NULL and PeaksLAp<>NULL
#print(Peakscommon)
} #if peak longer than max_peak
#print(Peakscommon)

} #for i 
#Peakscommon=Peakscommonold
#print(Peakscommon)
#write.table(Peakscommon,"lista_pików_330_10_tol0_0005_max_peak0_1_50.txt")

#print("Peakscommon:")
#print(Peakscommon)
##################The peaks of length <=10 removing
peak_length=10
Dlpeaks=NULL
for (yy in 1:nrow(Peakscommon)){
if (as.double(Peakscommon[yy,2])-as.double(Peakscommon[yy,1])>peak_length){
Dlpeaks=rbind(Dlpeaks,Peakscommon[yy,])}
}
Peakscommon=Dlpeaks

######################
Peakscommon=as.matrix(Peakscommon[which((Peakscommon[,1]>=longs[1])&(Peakscommon[,2]<=longs[2])),])
if (ncol(Peakscommon)==1) {Peakscommon=t(Peakscommon)}


if (length(Peakscommon)>=1){
##########3. The chromatograms integration over peaks for all peaks and all samples
# Labels - labels for samples in X
# The labels for the peaks are in the file lista_pików...txt
# TC - table with samples and peaks filled with integrals
TC=matrix(0,nrow(X),nrow(Peakscommon))
Y=t(X)
Y=abs(Y)  
x1=as.double(time[1:L])
for (i in 1:nrow(X)){
y1=as.double(Y[,i])
for (j in 1:nrow(Peakscommon)){
x=x1[Peakscommon[j,1]:Peakscommon[j,2]]
y=y1[Peakscommon[j,1]:Peakscommon[j,2]]
suma=sum(y)
TC[i,j]=suma
}
}

} #end if (!is.null(Peakscommon))
} #end if (length(Peakscommon)>=1)

#################
if (longs[3]<4.5){
text1=paste("pvclust_",wave,"_signal_",longs[3],"_",longs[4],"_gr",nrgr,"_peak_detection_part1_10_tol0_0005_max_peak0_1_50.pdf",sep="")
pdf(text1,width=10,height=6)
matplot(time[1:L],t(X),ylim=c(-0.06,0.06),xlim=c(2.5,4.5),type="l",main=paste("Dane po cow z pikami",wave," - part",sep=""))
if (!is.null(Peakscommon)) {abline(v=Peakscommon[,3])}
if (!is.null(Peakscommon)) {abline(v=Peakscommon[,4],col="red")}
dev.off()
}
if ((longs[3]>=4.5)&(longs[3]<=6.5)){
text2=paste("pvclust_",wave,"_signal_",longs[3],"_",longs[4],"_gr",nrgr,"_peak_detection_part2_10_tol0_0005_max_peak0_1_50.pdf",sep="")
pdf(text2,width=10,height=6)
matplot(time[1:L],t(X),ylim=c(-0.06,0.06),xlim=c(4.5,6.5),type="l",main=paste("Dane po cow z pikami",wave," - part",sep=""))
if (!is.null(Peakscommon)) {abline(v=Peakscommon[,3])}
if (!is.null(Peakscommon)) {abline(v=Peakscommon[,4],col="red")}
dev.off()
}
if ((longs[3]>=6.5)&(longs[3]<=8.5)){
text3=paste("pvclust_",wave,"_signal_",longs[3],"_",longs[4],"_gr",nrgr,"_peak_detection_part3_10_tol0_0005_max_peak0_1_50.pdf",sep="")
pdf(text3,width=10,height=6)
matplot(time[1:L],t(X),ylim=c(-0.06,0.06),xlim=c(6.5,8.5),type="l",main=paste("Dane po cow z pikami",wave," - part",sep=""))
if (!is.null(Peakscommon)) {abline(v=Peakscommon[,3])}
if (!is.null(Peakscommon)) {abline(v=Peakscommon[,4],col="red")}
dev.off()
}
if ((longs[3]>=8.5)&(longs[3]<=10.5)){
text4=paste("pvclust_",wave,"_signal_",longs[3],"_",longs[4],"_gr",nrgr,"_peak_detection_part4_10_tol0_0005_max_peak0_1_50.pdf",sep="")
pdf(text4,width=10,height=6)
matplot(time[1:L],t(X),ylim=c(-0.06,0.06),xlim=c(8.5,10.5),type="l",main=paste("Dane po cow z pikami",wave," - part",sep=""))
if (!is.null(Peakscommon)) {abline(v=Peakscommon[,3])}
if (!is.null(Peakscommon)) {abline(v=Peakscommon[,4],col="red")}
dev.off()
}
if (longs[3]>=10.5){
text5=paste("pvclust_",wave,"_signal_",longs[3],"_",longs[4],"_gr",nrgr,"_peak_detection_part5_10_tol0_0005_max_peak0_1_50.pdf",sep="")
pdf(text5,width=10,height=6)
matplot(time[1:L],t(X),ylim=c(-0.06,0.06),xlim=c(10.5,13),type="l",main=paste("Dane po cow z pikami",wave," - part",sep=""))
if (!is.null(Peakscommon)) {abline(v=Peakscommon[,3])}
if (!is.null(Peakscommon)) {abline(v=Peakscommon[,4],col="red")}
dev.off()
}

Peakscommon1<<-Peakscommon


if ((is.null(Peakscommon))|(length(Peakscommon)==0)){TC=0} else {TC}
}#end of function fun_peaks



Sspr=list(1:nrow(peaks))
peaksnew=NULL
c=NULL
table3=table[2:nrow(table),]
firstrow=table[1,]
countn=1
for (nrlongs in 1:nrow(peaks)){
if (as.double(peaks[nrlongs,7])>=0.1){
table2=table[2:nrow(table),1:8]
longs=as.double(peaks[nrlongs,3:6])
V=Xall[,longs[1]:longs[2]]
##
V=t(V)
colnames(V)=c(1:nrow(Xall))
cl <- pvclust(V,nboot=1000)
#plot(cl0, cex=0.85, cex.pv=0.7)
#ask.bak <- par()$ask
#par(ask=TRUE)
#pvrect(cl0,alpha=0.8)
cl1 <- pvpick(cl, alpha=0.85)
cl1
len=NULL
for (nrgr in 1:length(cl1$clusters)){
len=c(len,length(cl1$clusters[[nrgr]]))
}
print(len)
print(sort(len))

####Sample labels for each group importing to the file
tab01=NULL
tab02=NULL
tab00=table[2:nrow(table),]
rownames(tab00)=c(1:nrow(tab00))
for (oo in 1:length(cl1$clusters)){
tab01=cbind(tab00[cl1$clusters[[oo]],c(1,5:8)],rep(oo,length(cl1$clusters[[oo]])))
tab02=rbind(tab02,tab01)
}
write.csv(tab02,paste("pvclust_",wave,"_signal_",longs[3],"_",longs[4],"_nrlongs",nrlongs,"_labels.csv",sep=""))

####

##
TC1=list(1:length(cl1$clusters))
peaksnew=NULL
for (nrgr in 1:length(cl1$clusters)){
X=Xall[cl1$clusters[[nrgr]],]
TC1[[nrgr]]=fun_peaks(X,wave,time,L,longs)
pr=ncol(table2)
if ((!is.null(Peakscommon1))&(length(Peakscommon1)!=0)){table2=cbind(table2,matrix(0,nrow(table2),nrow(Peakscommon1)))
#if (ncol(TC1[[nrgr]])==1){TC1[[nrgr]]=t(TC1[[nrgr]])}
#table2 - TC in the right places, but without identification within tooo long peaks
#values are taken from table2, but plus info on the identification between groups within the peak
table2[as.double(cl1$clusters[[nrgr]]),(pr+1):(pr+nrow(Peakscommon1))]=TC1[[nrgr]]
POM=NULL
for (rr in 1:nrow(Peakscommon1)){
POM=rbind(POM,c(peaks[nrlongs,1:6],nrgr))
}
peaksnew=rbind(peaksnew,cbind(POM,Peakscommon1,matrix(0,nrow(Peakscommon1),1)))
} #end if !is.null(Peakscommon1)|(length(Peakscommon)!=0)
} #end nrgr


S=cbind(c(1:nrow(peaksnew)),peaksnew)
S1=NULL
for (r in 1:ncol(S)){
S1=cbind(S1,as.double(S[,r]))
}
S=S1
M=matrix(0,nrow(S),nrow(S))
for (w in 1:nrow(S)){
for (k in 1:nrow(S)){
if ((S[w,11]<=S[k,11])&(S[k,12]<=S[w,12])&(S[w,11]<=S[k,12])){
M[w,k]=(S[k,12]-S[k,11])/min(S[k,12]-S[k,11],S[w,12]-S[w,11])*100} else {
if ((S[k,11]<S[w,11])&(S[k,12]<=S[w,12])&(S[w,11]<=S[k,12])){
M[w,k]=(S[k,12]-S[w,11])/min(S[k,12]-S[k,11],S[w,12]-S[w,11])*100} else {
if ((S[w,11]<=S[k,11])&(S[w,12]<S[k,12])&(S[k,11]<=S[w,12])){
M[w,k]=(S[w,12]-S[k,11])/min(S[k,12]-S[k,11],S[w,12]-S[w,11])*100} else {
if ((S[k,11]<=S[w,11])&(S[w,12]<S[k,12])){
M[w,k]=(S[w,12]-S[w,11])/min(S[k,12]-S[k,11],S[w,12]-S[w,11])*100}
} #else
} #else
} #else
} #for
} #for

M=cbind(M,matrix(0,nrow(S),1))
for (w in 1:nrow(S)){
for (k in 1:nrow(S)){
if (M[w,k]>=80) {M[w,(nrow(S)+1)]=M[w,(nrow(S)+1)]+1}
}
}

for (w in 1:nrow(S)){
idg=which (M[w,]>=80)
idd=idg[which.max(M[which(M[w,(1:nrow(S))]>=80),(nrow(S)+1)])]
S[w,13]=idd
}


Sbiggerpeaks=function(S){
for (w in 1:nrow(S)){
if (S[S[w,13],13]!=S[w,13]){S[w,13]=S[S[w,13],13]}
}
S
}

Sv1=S[,13]
print(sort(unique(S[,13])))
S=Sbiggerpeaks(S)
Sv2=S[,13]
while (any(Sv1!=Sv2)) {Sv1=S[,13]
S=Sbiggerpeaks(S)
Sv2=S[,13]
print(sort(unique(S[,13])))}


su=sort(unique(S[,13])) #in S[,13] there are peaks with whom the peak is identified
#su - numbers of peaks selected to divide too long peak
#Sspr[[nrlongs]]=S
write.csv(S,paste("pvclust_",wave,"_signal_",longs[3],"_",longs[4],"_S_nrlongs",nrlongs,".csv",sep=""))

#Completed TC and a joint list of new peaks
Sr=matrix(0,length(su),6)
if (length(su)==1) {Sr=t(as.matrix(c(S[su,2:3],S[su,9:12])))} else {Sr=cbind(S[su,2:3],S[su,9:12])}
Sr[,1]=as.double(Sr[,1])*100+c(1:nrow(Sr))
pomfirstrow=NULL
for (pf in 1:nrow(Sr)){
pomfirstrow=c(pomfirstrow,paste(wave,"_",Sr[pf,5],"_",Sr[pf,6],sep=""))
}
print(nrlongs)
print(countn)
print(nrlongs+countn-2)
print(nrow(peakst))

if (nrlongs!=nrow(peaks)) {pomtp=peakst[(nrlongs+countn):nrow(peakst),]}
peakst=rbind(peakst[1:(nrlongs+countn-2),],Sr)
if (nrlongs!=nrow(peaks)) {peakst=rbind(peakst,pomtp)}  
print(peakst)

write.csv(Sr,paste("pvclust_",wave,"_peaks_Sr_nrlongs",nrlongs,".csv",sep=""), row.names=F)

#integrals are putting in the proper places after identification between groups within the peak
a=8
# if in TC1 elements are zero for chromatograms from outside the group, then columns can be added for the identified peak
# and then immediately there is the sum, if there is within group identification with the same peak

#table2i 
table2i=matrix(0,nrow(table2),length(su))

for (t2 in 1:nrow(S)) {table2i[,which(su==S[t2,13])]=table2i[,which(su==S[t2,13])]+table2[,t2+a]}

theend=ncol(table3)

if (nrlongs!=nrow(peaks)) {pomt3=table3[,(nrlongs+countn+a):theend]}
table3=cbind(table3[,1:(nrlongs+countn-2+a)],table2i)
if (nrlongs!=nrow(peaks)) {table3=cbind(table3,pomt3)}

write.csv(table2i,paste("pvclust_",wave,"_table2i_nrlongs",nrlongs,".csv",sep=""), row.names=F)

#(nrow(peakst)-length(su)+1+a)
if (nrlongs!=nrow(peaks)) {pomf=firstrow[(nrlongs+countn+a):theend]}
firstrow=c(firstrow[1:(nrlongs+countn-2+a)],pomfirstrow)
if (nrlongs!=nrow(piki)) {firstrow=c(firstrow,pomf)}

#wchodzi np dla nrlongs 14, 19
countn=countn+length(su)-1
} #end if >=0.1min
} #end nrlongs
firstr=as.vector(data.frame(firstrow),mode="character")
write.csv(peakst,paste("pvclust_",wave,"_peaks_indent.csv",sep=""))
tabel3=rbind(firstr,table3)
write.csv(table3,paste("pvclust_",wave,"_zad15_dosw3_phenotypedata.csv",sep=""))
