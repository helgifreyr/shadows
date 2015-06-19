(* ::Package:: *)

 Off[General::spell1]
Remove["Global`*"];
Unprotect[In,Out];
Clear[In,Out];



SetDirectory["/home/h/skoli/pt/blafis/phi-4"]
datOld = ReadList["global-data-BH-w=0.950.dat"];
lung1=Length[datOld];
 (*         1      2         3       4       5      6       7          8*)
(* Table[{w[[i]],rh[[i]],Mass[[i]],J[[i]],TH[[i]],S[[i]],Jint[[i]],Mint[[i]]},{i,1,lung0}] *)
wold=Table[  datOld[[i]][[19]],{i,   lung1 ,1,-1 }];
rhold=Table[ datOld[[i]][[1]],{i, lung1 ,1,-1 }];
Mass=Table[  datOld[[i]][[8]],{i, lung1 ,1,-1 }];
Jold=Table[  datOld[[i]][[20]],{i, lung1 ,1,-1 }]/2;
TH=Table[  datOld[[i]][[7]],{i,lung1 ,1,-1 }];
AH=Table[  datOld[[i]][[9]],{i,lung1 ,1,-1 }];
Q=Table[  datOld[[i]][[6]],{i,lung1 ,1,-1 }];


q=Q/Jold
aH=AH/(16 \[Pi] Mass^2);
tH=8 \[Pi] TH Mass;
j=-(Jold/ Mass^2);


dat = ReadList["global-data-multipoles-BH-w=0.950.dat"];

lung0=Length[dat];
 
 (*             1   2   3        4   5      6     7   8    9     10   11   12   13  14  15 
   Table[{  wf, ct,errFx0,errFxx0,ctF2,errFxx2,a5,ctF1,errFx1,b3F2,b3F1,b3F0,cfi,wt,rh}] 
  *) 

  w=Table[  dat[[i]][[1]],{i,1,lung0}];
ct= Table[  dat[[i]][[2]],{i,1,lung0}] ; 
errFx0= Table[  dat[[i]][[3]],{i,1,lung0}] ; 
errFxx0= Table[  dat[[i]][[4]],{i,1,lung0}] ; 
ctF2=  Table[  dat[[i]][[5]],{i,1,lung0}]; 
errFxx2=   Table[  dat[[i]][[6]],{i,1,lung0}];  
a5=  Table[  dat[[i]][[7]],{i,1,lung0}]; 
ctF1=  Table[  dat[[i]][[8]],{i,1,lung0}]; 
errFx1= Table[  dat[[i]][[9]],{i,1,lung0}]; 
b3F2= Table[  dat[[i]][[10]],{i,1,lung0}]; 
b3F1= Table[  dat[[i]][[11]],{i,1,lung0}]; 
b3F0= Table[  dat[[i]][[12]],{i,1,lung0}]; 
cfi= Table[  dat[[i]][[13]],{i,1,lung0}]; 
wt= Table[  dat[[i]][[14]],{i,1,lung0}]; 
rh= Table[  dat[[i]][[15]],{i,1,lung0}]; 


Qt=2 b3F0+a5 (2 ct-rh)+1/12 ct (-2 ct^2+3 ct rh+3 rh^2);
M=1/2 (-2 ct+rh);
J=1/2 cfi;
a=J/M^2;
rat=-((Qt M)/J^2);
jnew=J/M^2;



(*
 Finding equipotential surfaces for q

(J/M^2 , rat )  for given q;

*)
 
 (*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  I plot the full set of data and interpolate it 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
 

ini=1 ;
cut=0;
 
asa1=Table[{ q[[i]],j[[i]] },{i,ini,lung0-cut }];
asa2=Table[{ q[[i]],rat[[i]] },{i,ini,lung0-cut }];  


 Print["Plot (q,j)"]; 
ListPlot[asa1, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]

 Print["Plot (q,rat)"]; 
ListPlot[asa2, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


 (* I interpolate the existing data 
compute (J/M^2 , rat ) for a given value of M *)
 
ij=Interpolation[asa1,InterpolationOrder->3];
irat=Interpolation[asa2,InterpolationOrder->3]; 
 Print["total no of points = ",lung0];

q1=Table[{ q[[i]] },{i,1,lung0-cut }]; 
qmin=Min[q1];
qmax=Max[q1];

Print["max q=",qmax];
Print["min q=",qmin];


SetDirectory["multipoles"] ;
qs = Join[Table[i*0.05,{i,1,17}],{0.90,0.91,0.92,0.93,0.94,0.95,0.99,0.998}];
Do[
valq = qs[[i]];
If[qmin<valq<qmax,
dat = TableForm[{w[[1]],valq,ij[valq],irat[valq]}];
OutputForm[dat] >>> "q="<>ToString[valq]<>".txt";
]
,{i,1,Length[qs]}]
SetDirectory[".."];


(* I interpolate the existing data *)

asa=Table[{rat[[i]],M[[i]] },{i,1,lung0}];
ListPlot[asa, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]
iM=Interpolation[asa,InterpolationOrder->3];
qmin=Min[Abs[rat]];
qmax=Max[Abs[rat]];

Print["max q=",qmax];
Print["min q=",qmin];




v1=-25;
iM[v1]


Print["Grafic a-rat"];
s30=Table[{-a[[i]],-rat[[i]] },{i,1,lung0}];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


Print["Grafic rh-rat"];
s30=Table[{rh[[i]],-rat[[i]] },{i,1,lung0}];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


rh


Print["Grafic rh-ct"];
s30=Table[{rh[[i]],-ct[[i]] },{i,1,lung0}];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


Print["Grafic rh-cfi"];
s30=Table[{rh[[i]],-cfi[[i]] },{i,1,lung0}];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


Print["Grafic rh-a5"];
s30=Table[{rh[[i]],a5[[i]] },{i,1,lung0}];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


Print["Grafic rh-b3F0"];
s30=Table[{rh[[i]],b3F0[[i]] },{i,1, lung0 }];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


Print["Grafic rh-wt"];
s30=Table[{rh[[i]],wt[[i]] },{i,1,lung0 }];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]






(* now, study various errors *)


Print["Grafic rh-errFx0"];
s30=Table[{rh[[i]],errFx0[[i]] },{i,1,lung0}];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


Print["Grafic rh-errFxx0"];
s30=Table[{rh[[i]],errFxx0[[i]] },{i,1,lung0 }];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


Print["Grafic rh-errFxx2"];
s30=Table[{rh[[i]],errFxx2[[i]] },{i,1,lung0 }];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


(*
Print["Grafic rh-b3F1"];
s30=Table[{rh[[i]],b3F1[[i]] },{i,1,lung0 }];
gjaH=ListPlot[s30, PlotRange->All,Axes\[Rule]None, Frame\[Rule]True,Prolog\[Rule]AbsolutePointSize[4]];
*)
(*

Print["Grafic rh-b3F2"];
s30=Table[{rh[[i]],b3F2[[i]] },{i,1,lung0 }];
gjaH=ListPlot[s30, PlotRange->All,Axes\[Rule]None, Frame\[Rule]True,Prolog\[Rule]AbsolutePointSize[4]];*)


Print["Grafic rh-b3F2/b3F0"];
s30=Table[{rh[[i]],b3F2[[i]]/b3F0[[i]] },{i,1,lung0 }];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


Print["Grafic rh-b3F1/b3F0"];
s30=Table[{rh[[i]],b3F1[[i]]/b3F0[[i]] },{i,4,lung0 }];
gjaH=ListPlot[s30, PlotRange->All,Axes->None, Frame->True,Prolog->AbsolutePointSize[4]]


rh


(* relations  in the Kerr limit ;
wt=2/15 ct Sqrt[ct (ct-rh)] (2 ct^2-3 ct rh+rh^2);
b1=1/6 ct^2 (2 ct-3 rh);
a5=1/4 ct (ct-rh);
cfi=Sqrt[ct (ct-rh)] (-2 ct+rh);
*)


wtKerr=2/15 ct Sqrt[ct (ct-rh)] (2 ct^2-3 ct rh+rh^2);
b1Kerr=1/6 ct^2 (2 ct-3 rh);
a5Kerr=1/4 ct (ct-rh);
cfiKerr=Sqrt[ct (ct-rh)] (-2 ct+rh);


wt/wtKerr


b3F0/b1Kerr


a5/a5Kerr


cfi/cfiKerr
