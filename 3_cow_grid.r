#COW - correlation optimized warping at the grid of selected points


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
variety_ref=read.table("id_line_ref.txt") #the reference variety number loading
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

setwd(track_data_general)

#Table of observation import
LAB=read.csv("tabela_obserwacji_zad15_doœw3.csv")

obk=unique(LAB[,2])
obk[variety_ref] #reference variety

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

##########condition related to the loop id_ch, because in *.sh there is a condition - up to 24,
#because that is the maximum number of chromatograms for a line or variety, but of course there are lines in which there are fewer chromatograms
#if (id_ch<=ilm){


#matrices will include only selected variety
setwd(track_data_wave) 
matrices <- lapply(filenames[nrtab], read.table,header=T,fill=T)
setwd(track_results) 

#if (id_ch==1){
text00=paste("labels_samp_",id_wave,"_",id_variety,".txt",sep="")
write.table(Labels,text00)
#}

########Data normalization
#All chromatograms are dividing by the mass of each chromatogram
masa=NULL
for (mi in 1:ilm){
masa=c(masa,matrices[[mi]][1,3])
}

for (mi in 1:ilm){
matrices[[mi]][2402:15603,2]=matrices[[mi]][2402:15603,2]/masa[mi]
}

########Baseline removal based on differentiation:
for (mi in 1:ilm){
matrices[[mi]][2402:15602,2]=diff(matrices[[mi]][2402:15603,2])
}

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

########Function cow_opt_grid for the pair of chromatograms
###########################

fun_cow_opt_grid=function(L,time,ref,sec,lsec,mv,tv,id_wave){
###########################
#############!!!!!!!!New ref after warping is after differentation
d1=ref
d2=sec

############################
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
#ptl=Sys.time()
#print(ptl)
if (m>(t+3)) {
N=Lp/m

############################
#COW
T=d1
P=d2

############################################################################
time_o=c(1:L)
delta=0
LT=m*N
fsum=0
F=matrix(-Inf,(N+1),(LT+1))
U=F
F[N+1,N*m+1]=0
	
for (i in N:1) {
#cat("i:")
#print(i)
#ptl=Sys.time()
#print(ptl)
s1=max((i-1)*(m+delta-t),LT-(N-i+1)*(m+delta+t))+1

#cat("xstart:")
#print(s1)
e1=min((i-1)*(m+delta+t),LT-(N-i+1)*(m+delta-t))+1
#cat("xend:")
#print(e1)
for (x in s1:e1) { for (u in (delta-t):(delta+t)){
#dla róznych u countymy rózne F(i,x) i wybieramy max
#cat("x:")
#print(x)
##############
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
}
}
}
text0=as.character(lsec)
if (nchar(text0)==2){
text1=paste("grid",id_variety,"_",id_wave,"_",text0,sep="")} else {
text1=paste("grid",id_variety,"_",id_wave,"_0",text0,sep="")
}
text2=paste(text1,".txt",sep="")
write.table(R,text2)

print("The correct ending")

ptl=Sys.time()
print(ptl)
} #the end of funtion fun_cow_opt_grid


#Readout how much files are complited and filling the gaps:
filenamescheck <- dir(pattern=paste("grid",id_variety,"_",as.character(id_wave),sep=""),full.names=F)  
wyk=as.real(substr(filenamescheck,13,14))


########cow_opt_grid for whole variety
for (bb in 1:ilm){
if (!(any(wyk==bb))){ #if there is no file for the chromatogram, it it done
if (bb!=ind_ref){
sec=as.vector(matrices[[bb]][2402:15602,2])
lsec=bb
#function cow_opt for the pair of chromatograms
fun_cow_opt_grid(L,time,ref,sec,lsec,mv,tv,id_wave)
}
}
}

text0=as.character(ind_ref)
if (nchar(text0)==2){
text1=paste("grid",id_variety,"_",id_wave,"_",text0,sep="")} else {
text1=paste("grid",id_variety,"_",id_wave,"_0",text0,sep="")
}
text2=paste(text1,".txt",sep="")
write.table("0",text2)
#}
#do if id_chr<=24


