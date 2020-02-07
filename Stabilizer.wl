(* ::Package:: *)

BeginPackage["Stabilizer`"]


commE::usage = "commE[\[Sigma][i], \[Sigma][j]] checks if a pair of Pauli matrices \[Sigma][i], \[Sigma][j] commutes."
commE[\[Sigma][i_],\[Sigma][j_]]:=1/;(i==j || i==0 || j==0);
commE[\[Sigma][i_],\[Sigma][j_]]:=-1/;(i!=j && i!=0 && j!=0);


isComm::usage = "isComm[a_List, b_List] checks if a pair of tensor products of Pauli matrices commutes. a_List and b_List must be in the form of {\[Sigma][i1], \[Sigma][i2], ..., \[Sigma][iL]}."
isComm[a_List, b_List]:=Times@@MapThread[commE,{a,b}]


multiPauli::usage = "multiPauli[\[Sigma][i], \[Sigma][j]] returns a product of two Pauli matrices, \[Sigma][i], \[Sigma][j], up to multiplicative constant."
multiPauli[\[Sigma][i_],\[Sigma][j_]]:=(\[Sigma][i]\[Sigma][j]/.{\[Sigma][0]->1})/;(i*j==0&&i+j!=0);
multiPauli[\[Sigma][i_],\[Sigma][j_]]:=\[Sigma][0]/;(i==0&&j==0);
multiPauli[\[Sigma][i_],\[Sigma][j_]]:=\[Sigma][0]/;(i==j);
multiPauli[\[Sigma][i_],\[Sigma][j_]]:=\[Sigma]@@Complement[Range[3],{i,j}]/;(i!=j && i!=0 && j!=0);


Mat2Pauli::argerr="Valid input must be a {2L, 2L+1} dimensional matrix.";
Mat2Pauli[mat_, L_]:=Module[{sign,stabilizer,destabilizer},
If[Dimensions[mat] == {2L, 2L+1},
sign=((-1)^#&/@mat[[All,-1]]); stabilizer=Transpose/@MapThread[List,{mat[[All,;;L]],mat[[All,L+1;;-2]]}]/.{{0,0}->\[Sigma][0],{1,0}->\[Sigma][1],{1,1}->\[Sigma][2],{0,1}->\[Sigma][3]};
destabilizer=TensorProduct@@@stabilizer[[;;L]]*sign[[;;L]]; stabilizer=TensorProduct@@@stabilizer[[L+1;;]]*sign[[L+1;;]];
Return[{destabilizer, stabilizer}],
Message[Mat2Pauli::argerr]]];


CNOT::argerr="The matrix input must be {2L, 2L+1} dimensional. a \[LessEqual] L. b \[LessEqual] L.";
CNOT[mat_, L_, a_, b_]:=If[Dimensions[mat] == {2L, 2L+1} && a<=L && b<=L,
Module[{res=mat},
res[[All,-1]]=Mod[Mod[mat[[All,b]]+mat[[All, L+a]]+1,2] mat[[All,a]]mat[[All,L+b]] + mat[[All,-1]],2];
res[[All,b]]=Mod[mat[[All,b]]+mat[[All,a]],2];
res[[All,L+a]]= Mod[mat[[All,L+b]]+mat[[All,L+a]],2];
Return[res]],
Message[CNOT::argerr]];


PHASE::argerr="The matrix input must be {2L, 2L+1} dimensional. a \[LessEqual] L.";
PHASE[mat_, L_, a_]:=If[Dimensions[mat] == {2L, 2L+1} && a<=L,
Module[{res=mat},
res[[All,-1]]=Mod[mat[[All,a]]mat[[All,L+a]]+mat[[All,-1]],2];
res[[All,L+a]]= Mod[mat[[All,L+a]]+mat[[All,L]],2];
Return[res]],
Message[PHASE::argerr]];


HAM::argerr="The matrix input must be {2L, 2L+1} dimensional. a \[LessEqual] L.";
HAM[mat_, L_, a_]:=If[Dimensions[mat] == {2L, 2L+1} && a<=L,
Module[{res=mat},
res[[All,-1]]=Mod[mat[[All,a]]mat[[All,L+a]]+mat[[All,-1]],2];
res[[All,L+a]]= mat[[All,a]]; res[[All,a]]= mat[[All,a+L]];
Return[res]],
Message[HAM::argerr]];


EndPackage[]
