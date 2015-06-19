(* ::Package:: *)

 Off[General::spell1]
Remove["Global`*"];
Unprotect[In,Out];
Clear[In,Out];


dat = ReadList["m=1.0-c2=350.0-w=0.850.dat"];


lung0=Length[dat]

 (*       1   2   3  4  5   6     7  8   9  10   11   12  13   14  15  16 17 18 19  20 Table[{rh,alfa,c1,c2,c3,Jint,THc,Mc,AHc,err1,minF0,f0H,f1H,f2H,ZH,Mint,Le,Lp,w,constJINF,
  21     22   23    24     25   26   27    28   29   30
maxF0,minF1,maxF1,minF2,maxF2,minW,maxW,minZ,maxZ,errSmarr}]   *) 

rh=Table[  dat[[i]][[1]],{i,1,lung0}];
alfa=Table[  dat[[i]][[2]],{i,1,lung0}];
c1= Table[  dat[[i]][[3]],{i,1,lung0}] ;(* coefficient Z^6 in the potential *)
c2= Table[  dat[[i]][[4]],{i,1,lung0}] ;(* coefficient Z^4 in the potential *)
c3= Table[  dat[[i]][[5]],{i,1,lung0}] ;(* coefficient Z^2 in the potential -- mass term *)
Jint=Table[  dat[[i]][[6]],{i,1,lung0}];
THc=  Table[  dat[[i]][[7]],{i,1,lung0}]; 
Mc=  Table[  dat[[i]][[8]],{i,1,lung0}]; 
AHc= Table[  dat[[i]][[9]],{i,1,lung0}]; 
 F0h= Table[  dat[[i]][[12]],{i,1,lung0}]; 
 F1h= Table[  dat[[i]][[13]],{i,1,lung0}]; 
 F2h= Table[  dat[[i]][[14]],{i,1,lung0}]; 
Zh= Table[  dat[[i]][[15]],{i,1,lung0}]; 
Mint= Table[  dat[[i]][[16]],{i,1,lung0}];
Le= Table[  dat[[i]][[17]],{i,1,lung0}];
Lp= Table[  dat[[i]][[18]],{i,1,lung0}];
w= Table[  dat[[i]][[19]],{i,1,lung0}];
J=1/2 Table[  dat[[i]][[20]],{i,1,lung0}];
minF0= Table[  dat[[i]][[11]],{i,1,lung0}]; 
maxF0= Table[  dat[[i]][[21]],{i,1,lung0}]; 
minF1= Table[  dat[[i]][[22]],{i,1,lung0}]; 
maxF1= Table[  dat[[i]][[23]],{i,1,lung0}]; 
minF2= Table[  dat[[i]][[24]],{i,1,lung0}]; 
maxF2= Table[  dat[[i]][[25]],{i,1,lung0}]; 
minW= Table[  dat[[i]][[26]],{i,1,lung0}]; 
maxW= Table[  dat[[i]][[27]],{i,1,lung0}]; 
minZ= Table[  dat[[i]][[28]],{i,1,lung0}]; 
maxZ= Table[  dat[[i]][[29]],{i,1,lung0}]; 
errSmarr= Table[  dat[[i]][[30]],{i,1,lung0}]; 


(* I define the Schw BH quantities & global data *)
AH0=4 Pi rh^2;
AH=AH0 AHc;

TH0=1/(4 Pi rh);
TH=TH0 THc;

M0=rh/2  ;
Mass=M0+Mc;
S=1/4 AH ;


Max[rh]


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 First I plot the full set of data;
then I split it into two parts -- the two branches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

Q=Jint;

q=Q/J;
nr1=Ordering[q,1][[1]];


asa=Table[{q[[i]],Mass[[i]] },{i,1,lung0}]; 
asa1=Table[{q[[i]],-J[[i]] },{i,1,lung0}]; 
asa2=Table[{q[[i]],AH[[i]] },{i,1,lung0}]; 
asa3=Table[{q[[i]],TH[[i]] },{i,1,lung0}]; 

qmin=Min[Abs[q]];
qmax=Max[Abs[q]];

Print["max q=",qmax];
Print["min q=",qmin];



(* I interpolate the existing data 
compute Mass for a given value of Q/J *)

asa=Table[{q[[i]],Mass[[i]] },{i,1,nr1}]; 
asa1=Table[{q[[i]],-J[[i]] },{i,1,nr1}]; 
asa2=Table[{q[[i]],AH[[i]] },{i,1,nr1}]; 
asa3=Table[{q[[i]],TH[[i]] },{i,1,nr1}]; 

iM=Interpolation[asa,InterpolationOrder->3];
iJ=Interpolation[asa1,InterpolationOrder->3];
iAH=Interpolation[asa2,InterpolationOrder->3];
iTH=Interpolation[asa3,InterpolationOrder->3];

q1=Table[{q[[i]]  },{i, 1,lung0}]; 
qmin=Min[Abs[q1]];
qmax=Max[Abs[q1]];

Print["max q=",qmax];
Print["min q=",qmin];





valq=0.985;
  Do[
dat1=TableForm[{{w[[1]],valq,iM[valq],iJ[valq],iAH[valq],iTH[valq]}}];
Print[dat1];
 OutputForm[dat1] >>> "q=0,985.txt";
 ]; 


valq=0.97; 
 Do[
dat1=TableForm[{{w[[1]],valq,iM[valq],iJ[valq],iAH[valq],iTH[valq]}}];
Print[dat1];
 OutputForm[dat1] >>> "q=0,97.txt";
 ]; 


(* I interpolate the existing data 
compute Mass for a given value of Q/J *)



asa=Table[{q[[i]],Mass[[i]] },{i,nr1+1,lung0}]; 
asa1=Table[{q[[i]],-J[[i]] },{i,nr1+1,lung0}]; 
asa2=Table[{q[[i]],AH[[i]] },{i,nr1+1,lung0}]; 
asa3=Table[{q[[i]],TH[[i]] },{i,nr1+1,lung0}]; 

iM=Interpolation[asa,InterpolationOrder->3];
iJ=Interpolation[asa1,InterpolationOrder->3];
iAH=Interpolation[asa2,InterpolationOrder->3];
iTH=Interpolation[asa3,InterpolationOrder->3];

qmin=Min[Abs[q]];
qmax=Max[Abs[q]];

Print["max q=",qmax];
Print["min q=",qmin];



valq=0.985;
Do[
dat1=TableForm[{{w[[1]],valq,iM[valq],iJ[valq],iAH[valq],iTH[valq]}}];
Print[dat1];
 OutputForm[dat1] >>> "q=0,985.txt";
 ]; 


valq=0.97;
Do[
dat1=TableForm[{{w[[1]],valq,iM[valq],iJ[valq],iAH[valq],iTH[valq]}}];
Print[dat1];
 OutputForm[dat1] >>> "q=0,97.txt";
 ]; 




