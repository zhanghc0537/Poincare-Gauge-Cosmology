(* ::Package:: *)

(************************ 0. Info and copyright ***********************)


PGC122`$Version={"1.2.2",{2020,6,30}}


(* PGC, Poincare Gauge Cosmology *)

(* Copyright (C) 2018-2020 Hongchao Zhang *)

(* Dalian University of Technology & Penn State University *)

(* This program is based on the xAct series *)

(* Version 1.2.2: Systematic approach, based on v0.5.3, accomplished *)


(************************ 1. Begin package ***********************)


BeginPackage["PGC122`",{"xAct`xCore`","xAct`xPerm`","xAct`xTensor`","xAct`xCoba`","xAct`xPert`"}]


(* Parallelization setting up *)
ParallelNeeds["xAct`xCore`"];
ParallelNeeds["xAct`xPerm`"];
ParallelNeeds["xAct`xTensor`"];
ParallelNeeds["xAct`xCoba`"];
ParallelNeeds["xAct`xPert`"];


PGC122`$barslength=60;
PGC122`bars=StringJoin[Table["-",{PGC122`$barslength}]];


Print[PGC122`bars];
Print[PGC122`bars];
Print["Package PGC`  version ",PGC122`$Version[[1]],", ",PGC122`$Version[[2]]];
Print["PGC, Poincare Gauge Cosmology"];
Print["Copyright (C) 2018-2020 Hongchao Zhang"];
Print["Dalian University of Technology & Penn State University"];
Print["This program is based on the xAct series"];
Print["Version 1.2.2: Systematic approach, based on v0.5.3, accomplished"];
Print[PGC122`bars];


(*usage*)


(*Begin["`Private`"]*)


(************************ 2. General Statements ***********************)


(*Block[{Print},
<<xAct`xTensor`;
ParallelNeeds["xAct`xTensor`"];
(*Needs["xAct`xPrint`"]*)
<<xAct`xCoba`;
<<xAct`xPert`;
];*)


$DefInfoQ=False;
$CVVerbose=False;
$PrePrint=ScreenDollarIndices;


(* Adding the working directory *)
$PackageDirectory=ParentDirectory[DirectoryName[FindFile["PGC122`"]]];
If[!ContainsAny[$Path,{$PackageDirectory}],PrependTo[$Path,$PackageDirectory]];
Print["The current working directory is: ",$PackageDirectory];
Print[PGC122`bars];


(* my options: Torsion or Contorsion *)
$WhichVar=1;
(* Pre evaluating cdcdTorsionCD or not: 1 is yes *)
$PreEvaluatecdcdTorsionCD=1;


(* Define a manifold M *)
DefManifold[M,4,{\[Mu],\[Nu],\[Rho],\[Sigma],\[Tau],\[Lambda],\[Upsilon],\[Xi],\[Zeta],\[Alpha],\[Beta],\[Gamma],i,j,k,l,m,n,a,b,c,d,e},PrintAs->"\!\(\*SuperscriptBox[\(M\), \(4\)]\)"];


(* Define a metric g and its perturbation *)
DefMetric[-1,metric[-\[Mu],-\[Nu]],cd,SymbolOfCovD->{";","D"},WeightedWithBasis->AIndex,PrintAs->"g"];
DefMetricPerturbation[metric,metpert,\[Epsilon]];
PrintAs[metpert]^="\[Delta]g";


Unprotect[IndexForm];
IndexForm[LI[x_]]:=ColorString[ToString[x],RGBColor[0,0,1]];
Protect[IndexForm];


(* Define a metric-affine connection compatible with the general covariant derivative CD and metric g *)
DefCovD[CD[-\[Mu]],Torsion->True,SymbolOfCovD->{"|","\[Del]"},FromMetric->metric];


(*Define the perturbations of metric-affine connection \[CapitalGamma], torsion T, and contorsion K respectively *)
DefTensorPerturbation[PertConnection[LI[order],\[Sigma],-\[Mu],-\[Nu]],ChristoffelCD[\[Sigma],-\[Mu],-\[Nu]],M];
PrintAs[PertConnection]^="\[Delta]\[CapitalGamma]";
DefTensorPerturbation[PertTorsion[LI[order],\[Sigma],-\[Mu],-\[Nu]],TorsionCD[\[Sigma],-\[Mu],-\[Nu]],M];
PrintAs[PertTorsion]^="\[Delta]T";
DefTensor[Contorsion[-\[Sigma],-\[Mu],-\[Nu]],M,Antisymmetric[{-\[Sigma],-\[Nu]}],PrintAs->"K"];
DefTensorPerturbation[PertContorsion[LI[order],-\[Sigma],-\[Mu],-\[Nu]],Contorsion[-\[Sigma],-\[Mu],-\[Nu]],M];
PrintAs[PertContorsion]^="\[Delta]K";


(* Rules among \[CapitalGamma], T and K *)
rule\[CapitalGamma]=ChristoffelCD[\[Sigma]_,\[Mu]_,\[Nu]_]->Christoffelcd[\[Sigma],\[Mu],\[Nu]]+Contorsion[\[Sigma],\[Mu],\[Nu]];
ruleK=Contorsion[\[Sigma]_,\[Mu]_,\[Nu]_]->1/2 (TorsionCD[\[Sigma],\[Mu],\[Nu]]+TorsionCD[\[Mu],\[Sigma],\[Nu]]+TorsionCD[\[Nu],\[Sigma],\[Mu]]);
(*rule=Which[$WhichVar==1,rule\[CapitalGamma]/.ruleK,$WhichVar==2,rule\[CapitalGamma]];*)


metorder=metpert[LI[order_],__]:>0/;order>1;


DefTensor[VTorsion[\[Sigma]],M,PrintAs->"V"];
DefTensor[ATorsion[\[Sigma]],M,PrintAs->"A"];
DefTensor[TTorsion[-\[Xi],-\[Mu],-\[Nu]],M,Symmetric[{-\[Xi],-\[Mu]}],PrintAs->"t"];
DefTensorPerturbation[PertVTorsion[LI[order],\[Sigma]],VTorsion[\[Sigma]],M,PrintAs->"\[Delta]V"];
DefTensorPerturbation[PertATorsion[LI[order],\[Sigma]],ATorsion[\[Sigma]],M,PrintAs->"\[Delta]A"];
DefTensorPerturbation[PertTTorsion[LI[order],-\[Sigma],-\[Mu],-\[Nu]],TTorsion[-\[Sigma],-\[Mu],-\[Nu]],M,Symmetric[{-\[Sigma],-\[Mu]}],PrintAs->"\[Delta]t"];
TTorsion/:TTorsion[\[Sigma]_,ChangeIndex@\[Sigma]_,\[Mu]_]:=0;
TTorsion/:TTorsion[ChangeIndex@\[Sigma]_,\[Sigma]_,\[Mu]_]:=0;
TTorsion/:TTorsion[\[Sigma]_,\[Mu]_,ChangeIndex@\[Sigma]_]:=0;
TTorsion/:TTorsion[ChangeIndex@\[Sigma]_,\[Mu]_,\[Sigma]_]:=0;
TTorsion/:TTorsion[\[Mu]_,\[Sigma]_,ChangeIndex@\[Sigma]_]:=0;
TTorsion/:TTorsion[\[Mu]_,ChangeIndex@\[Sigma]_,\[Sigma]_]:=0;
DefTensor[ACurvature[-\[Sigma],-\[Tau],-\[Mu],-\[Nu]],M,PrintAs->"A"];
DefTensor[BCurvature[-\[Sigma],-\[Tau],-\[Mu],-\[Nu]],M,PrintAs->"B"];
DefTensor[WCurvature[-\[Sigma],-\[Tau],-\[Mu],-\[Nu]],M,PrintAs->"W"];
DefTensor[CCurvature[-\[Sigma],-\[Tau],-\[Mu],-\[Nu]],M,PrintAs->"C"];
DefTensor[ICurvature[-\[Sigma],-\[Tau]],M,PrintAs->"I"];
DefTensor[ECurvature[-\[Sigma],-\[Tau]],M,PrintAs->"E"];


rulevt=VTorsion[\[Nu]_]:>1/3 Module[{\[Sigma]},TorsionCD[\[Sigma],\[Nu],ChangeIndex@\[Sigma]]];
ruleat=ATorsion[\[Nu]_]:>1/6 Module[{\[Xi],\[Mu],\[Zeta]},epsilonmetric[\[Nu],\[Xi],\[Mu],\[Zeta]]TorsionCD[ChangeIndex[\[Xi]],ChangeIndex[\[Mu]],ChangeIndex[\[Zeta]]]];
rulett=TTorsion[\[Mu]_,\[Nu]_,\[Sigma]_]:>1/2 (TorsionCD[\[Mu],\[Nu],\[Sigma]]+TorsionCD[\[Nu],\[Mu],\[Sigma]])-1/2 (metric[\[Sigma],\[Mu]]VTorsion[\[Nu]]+metric[\[Sigma],\[Nu]]VTorsion[\[Mu]])+metric[\[Mu],\[Nu]]VTorsion[\[Sigma]];
ruleac=ACurvature[\[Sigma]_,\[Tau]_,\[Mu]_,\[Nu]_]:>1/6 (RiemannCD[\[Sigma],\[Tau],\[Mu],\[Nu]]+RiemannCD[\[Sigma],\[Mu],\[Nu],\[Tau]]+RiemannCD[\[Sigma],\[Nu],\[Tau],\[Mu]]+RiemannCD[\[Tau],\[Mu],\[Sigma],\[Nu]]+RiemannCD[\[Tau],\[Nu],\[Mu],\[Sigma]]+RiemannCD[\[Mu],\[Nu],\[Sigma],\[Tau]]);
rulebc=BCurvature[\[Sigma]_,\[Tau]_,\[Mu]_,\[Nu]_]:>1/4 (WeylCD[\[Sigma],\[Tau],\[Mu],\[Nu]]+WeylCD[\[Mu],\[Nu],\[Sigma],\[Tau]]-WeylCD[\[Sigma],\[Nu],\[Tau],\[Mu]]-WeylCD[\[Tau],\[Mu],\[Sigma],\[Nu]]);
rulecc=CCurvature[\[Sigma]_,\[Tau]_,\[Mu]_,\[Nu]_]:>1/2 (WeylCD[\[Sigma],\[Tau],\[Mu],\[Nu]]-WeylCD[\[Mu],\[Nu],\[Sigma],\[Tau]]);
ruleic=ICurvature[\[Sigma]_,\[Tau]_]:>1/2 (RicciCD[\[Sigma],\[Tau]]+RicciCD[\[Tau],\[Sigma]])-1/4 metric[\[Sigma],\[Tau]]RicciScalarCD[];
ruleec=ECurvature[\[Sigma]_,\[Tau]_]:>1/2 (RicciCD[\[Sigma],\[Tau]]-RicciCD[\[Tau],\[Sigma]]);
ruleIrrToTC={rulevt,ruleat,rulett,ruleac,rulebc,rulecc,ruleic,ruleec};
ruleIrrDecomp=TorsionCD[\[Sigma]_,\[Mu]_,\[Nu]_]:>Module[{\[Zeta]},2/3 (TTorsion[\[Sigma],\[Mu],\[Nu]]-TTorsion[\[Sigma],\[Nu],\[Mu]])-(metric[\[Sigma],\[Mu]]VTorsion[\[Nu]]-metric[\[Sigma],\[Nu]]VTorsion[\[Mu]])+epsilonmetric[\[Sigma],\[Mu],\[Nu],\[Zeta]]ATorsion[ChangeIndex@\[Zeta]]];


Clear[b0,a1,a2,a3,b1,b2,b3,b4,b5,b6];
Clear[c1];
(*DefConstantSymbol[{b,a1,a2,a3,b1,b2,b3,b4,b5,b6}];
coef={b,a1,a2,a3,b1,b2,b3,b4,b5,b6};*)
(*DefConstantSymbol[\[Alpha]];
DefConstantSymbol[{a1,a2,a3}];
DefConstantSymbol[{b1,b2,b3,b4,b5,b6}];
DefConstantSymbol[c1];
DefScalarFunction[a];
DefScalarFunction[h];
DefScalarFunction[f];
DefScalarFunction[\[Phi]];
DefScalarFunction[B];
DefScalarFunction[\[Psi]];
DefScalarFunction[\[CapitalEpsilon]];
DefScalarFunction[\[Phi]2];
DefScalarFunction[B2];
DefScalarFunction[\[Psi]2];
DefScalarFunction[\[CapitalEpsilon]2];
DefScalarFunction[v];
DefScalarFunction[v2];
DefScalarFunction[\[CapitalXi]];
DefScalarFunction[\[CapitalTheta]];
DefScalarFunction[\[CapitalPhi]];
DefScalarFunction[\[CapitalPsi]];
DefScalarFunction[\[CapitalXi]2];
DefScalarFunction[\[CapitalTheta]2];
DefScalarFunction[\[CapitalPhi]2];
DefScalarFunction[\[CapitalPsi]2];
DefScalarFunction[\[Rho],PrintAs->"\!\(\*OverscriptBox[\(\[Rho]\), \(_\)]\)"];
DefScalarFunction[\[Rho]m,PrintAs->"\!\(\*SubscriptBox[OverscriptBox[\(\[Rho]\), \(_\)], \(m\)]\)"];
DefScalarFunction[\[Rho]r,PrintAs->"\!\(\*SubscriptBox[OverscriptBox[\(\[Rho]\), \(_\)], \(r\)]\)"];
DefScalarFunction[p,PrintAs->"\!\(\*OverscriptBox[\(p\), \(_\)]\)"];
DefConstantSymbol[w];
DefScalarFunction[H];
DefConstantSymbol[kappa,PrintAs->"\[Kappa]"];*)


DefBasis[tetrad,TangentM,{0,1,2,3},BasisColor->RGBColor[1,0,1]];
minkowski=DiagonalMatrix[{-1,1,1,1}];
AllComponentValues[Basis[{-\[Mu],-tetrad},{-\[Nu],-tetrad}],minkowski];

DefChart[chart,M,{0,1,2,3},{t[],x[],y[],z[]},ChartColor->RGBColor[1,0,0]];
AllComponentValues[metpert[LI[0],{-\[Mu],-chart},{-\[Nu],-chart}],minkowski];
MetricCompute[metric,chart,All,CVSimplify->Simplify];


DefTensor[EinsteinEqleft[\[Mu],\[Nu]],M,Symmetric[{\[Mu],\[Nu]}]];
DefTensor[CartanEqleft[\[Sigma],-\[Mu],-\[Nu]],M,Antisymmetric[{-\[Mu],-\[Nu]}]];
DefTensor[EinsteinEqleftK[\[Mu],\[Nu]],M,Symmetric[{\[Mu],\[Nu]}]];
DefTensor[CartanEqleftK[-\[Sigma],-\[Mu],-\[Nu]],M,Antisymmetric[{-\[Sigma],-\[Nu]}]];


DefConstantSymbol[{b0,a1,a2,a3,b1,b2,b3,b4,b5,b6}];
coef={b0,a1,a2,a3,b1,b2,b3,b4,b5,b6};
LagrangianT2=Module[{\[Sigma],\[Mu],\[Nu],\[Tau],\[Xi]},a1 TTorsion[\[Sigma],\[Mu],\[Nu]]TTorsion[ChangeIndex@\[Sigma],ChangeIndex@\[Mu],ChangeIndex@\[Nu]]+a2 VTorsion[\[Tau]]VTorsion[ChangeIndex@\[Tau]]+a3 ATorsion[\[Xi]]ATorsion[ChangeIndex@\[Xi]]];
LagrangianR2=Module[{\[Sigma],\[Tau],\[Mu],\[Nu]},b1 ACurvature[\[Sigma],\[Tau],\[Mu],\[Nu]]ACurvature[ChangeIndex@\[Sigma],ChangeIndex@\[Tau],ChangeIndex@\[Mu],ChangeIndex@\[Nu]]+b2 BCurvature[\[Sigma],\[Tau],\[Mu],\[Nu]]BCurvature[ChangeIndex@\[Sigma],ChangeIndex@\[Tau],ChangeIndex@\[Mu],ChangeIndex@\[Nu]]+b3 CCurvature[\[Sigma],\[Tau],\[Mu],\[Nu]]CCurvature[ChangeIndex@\[Sigma],ChangeIndex@\[Tau],ChangeIndex@\[Mu],ChangeIndex@\[Nu]]+b4 ECurvature[\[Mu],\[Nu]]ECurvature[ChangeIndex@\[Mu],ChangeIndex@\[Nu]]+b5 ICurvature[\[Mu],\[Nu]]ICurvature[ChangeIndex@\[Mu],ChangeIndex@\[Nu]]+b6 RicciScalarCD[]^2];
Lagrangian=b0 RicciScalarCD[]+LagrangianT2+LagrangianR2;


DefConstantSymbol[{A1,A2,A3,B1,B2,B3,B4,B5}];


DefTensor[\[Omega][\[Mu],\[Nu]],M,Symmetric[{\[Mu],\[Nu]}]];
DefTensor[\[Theta][\[Mu],\[Nu]],M,Symmetric[{\[Mu],\[Nu]}]];
\[Omega]/:\[Omega][\[Mu]_,\[Nu]_]\[Omega][ChangeIndex[\[Nu]_],\[Sigma]_]:=\[Omega][\[Mu],\[Sigma]];
\[Omega]/:\[Omega][\[Mu]_,\[Nu]_]\[Omega][\[Sigma]_,ChangeIndex[\[Nu]_]]:=\[Omega][\[Mu],\[Sigma]];
\[Theta]/:\[Theta][\[Mu]_,\[Nu]_]\[Theta][ChangeIndex[\[Nu]_],\[Sigma]_]:=\[Theta][\[Mu],\[Sigma]];
\[Theta]/:\[Theta][\[Mu]_,\[Nu]_]\[Theta][\[Sigma]_,ChangeIndex[\[Nu]_]]:=\[Theta][\[Mu],\[Sigma]];
\[Omega]/:\[Omega][\[Mu]_,\[Nu]_]\[Theta][ChangeIndex[\[Nu]_],\[Sigma]_]:=0;
\[Omega]/:\[Omega][\[Mu]_,\[Nu]_]\[Theta][\[Sigma]_,ChangeIndex[\[Nu]_]]:=0;
\[Theta]/:\[Theta][\[Mu]_,\[Nu]_]\[Omega][ChangeIndex[\[Nu]_],\[Sigma]_]:=0;
\[Theta]/:\[Theta][\[Mu]_,\[Nu]_]\[Omega][\[Sigma]_,ChangeIndex[\[Nu]_]]:=0;
\[Omega]/:\[Omega][\[Mu]_,ChangeIndex[\[Mu]_]]:=1;
\[Theta]/:\[Theta][\[Mu]_,ChangeIndex[\[Mu]_]]:=3;
DefTensor[kv[\[Mu]],M,PrintAs->"k"];
(*DefScalarFunction[kv2];*)
DefTensor[kv2[],M];
kv/:kv[\[Xi]_]kv[ChangeIndex[\[Xi]_]]:=kv2[];
DefTensor[k1[\[Mu]],M,PrintAs->"\!\(\*OverscriptBox[\(k\), \(\[Tilde]\)]\)"];
k1/:k1[\[Xi]_]k1[ChangeIndex[\[Xi]_]]:=1;
k1/:k1[\[Xi]_]k1[\[Zeta]_]:=\[Omega][\[Xi],\[Zeta]];
k1/:k1[\[Xi]_]\[Theta][ChangeIndex[\[Xi]_],\[Zeta]_]:=0;
k1/:k1[\[Xi]_]\[Theta][\[Zeta]_,ChangeIndex[\[Xi]_]]:=0;
k1/:k1[\[Xi]_]\[Omega][ChangeIndex[\[Xi]_],\[Zeta]_]:=k1[\[Zeta]];
k1/:k1[\[Xi]_]\[Omega][\[Zeta]_,ChangeIndex[\[Xi]_]]:=k1[\[Zeta]];
ruleProj={\[Omega][\[Mu]_,\[Nu]_]:>(kv[\[Mu]]kv[\[Nu]])/kv2[],\[Theta][\[Mu]_,\[Nu]_]:>metric[\[Mu],\[Nu]]-(kv[\[Mu]]kv[\[Nu]])/kv2[],k1[\[Mu]_]:>kv[\[Mu]]/Sqrt[kv2[]]};
(* definition for spin-0^- *)
DefTensor[PAA0m[i,j,k,c,a,b],M,Antisymmetric[{i,j,k}],PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(0\), \(A\), \(-\)]]\)"]
(* definition for spin-0^+ *)
DefTensor[Pat0p[i,j,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(a\), SubsuperscriptBox[\(0\), \(t\), \(+\)]]\)"];
DefTensor[Pats0p[i,j,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(a\), SubsuperscriptBox[\(0\), \(ts\), \(+\)]]\)"];DefTensor[PaAtV0p[i,j,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(aA\), SubsuperscriptBox[\(0\), \(tV\), \(+\)]]\)"];DefTensor[Past0p[i,j,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(a\), SubsuperscriptBox[\(0\), \(st\), \(+\)]]\)"];DefTensor[Pas0p[i,j,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(a\), SubsuperscriptBox[\(0\), \(s\), \(+\)]]\)"];DefTensor[PaAsV0p[i,j,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(aA\), SubsuperscriptBox[\(0\), \(sV\), \(+\)]]\)"];DefTensor[PAaVt0p[i,j,k,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(Aa\), SubsuperscriptBox[\(0\), \(Vt\), \(+\)]]\)"];DefTensor[PAaVs0p[i,j,k,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(Aa\), SubsuperscriptBox[\(0\), \(Vs\), \(+\)]]\)"];DefTensor[PAV0p[i,j,k,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(0\), \(V\), \(+\)]]\)"];
(* definition for spin-1^- *)
DefTensor[Pas1m[i,j,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(a\), SubsuperscriptBox[\(1\), \(s\), \(-\)]]\)"];
DefTensor[Pasm1m[i,j,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(a\), SubsuperscriptBox[\(1\), \(sm\), \(-\)]]\)"];
DefTensor[PaAsV1m[i,j,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(aA\), SubsuperscriptBox[\(1\), \(sV\), \(-\)]]\)"];
DefTensor[PaAst1m[i,j,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(aA\), SubsuperscriptBox[\(1\), \(st\), \(-\)]]\)"];
DefTensor[Pams1m[i,j,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(a\), SubsuperscriptBox[\(1\), \(ms\), \(-\)]]\)"];
DefTensor[Pam1m[i,j,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(a\), SubsuperscriptBox[\(1\), \(m\), \(-\)]]\)"];
DefTensor[PaAmV1m[i,j,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(aA\), SubsuperscriptBox[\(1\), \(mV\), \(-\)]]\)"];
DefTensor[PaAmt1m[i,j,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(aA\), SubsuperscriptBox[\(1\), \(mt\), \(-\)]]\)"];
DefTensor[PAaVs1m[i,j,k,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(Aa\), SubsuperscriptBox[\(1\), \(Vs\), \(-\)]]\)"];
DefTensor[PAaVm1m[i,j,k,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(Aa\), SubsuperscriptBox[\(1\), \(Vm\), \(-\)]]\)"];
DefTensor[PAV1m[i,j,k,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(1\), \(V\), \(-\)]]\)"];
DefTensor[PAVt1m[i,j,k,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(1\), \(Vt\), \(-\)]]\)"];
DefTensor[PAats1m[i,j,k,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(Aa\), SubsuperscriptBox[\(1\), \(ts\), \(-\)]]\)"];
DefTensor[PAatm1m[i,j,k,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(Aa\), SubsuperscriptBox[\(1\), \(tm\), \(-\)]]\)"];
DefTensor[PAtV1m[i,j,k,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(1\), \(tV\), \(-\)]]\)"];
DefTensor[PAt1m[i,j,k,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(1\), \(t\), \(-\)]]\)"];
(* definition for spin-1^+ *)
DefTensor[Pae1p[i,j,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(a\), SubsuperscriptBox[\(1\), \(e\), \(+\)]]\)"];
DefTensor[PaAeA1p[i,j,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(aA\), SubsuperscriptBox[\(1\), \(eA\), \(+\)]]\)"];
DefTensor[PaAet1p[i,j,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(aA\), SubsuperscriptBox[\(1\), \(et\), \(+\)]]\)"];
DefTensor[PAaAe1p[i,j,k,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(Aa\), SubsuperscriptBox[\(1\), \(Ae\), \(+\)]]\)"];
DefTensor[PAA1p[i,j,k,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(1\), \(A\), \(+\)]]\)"];
DefTensor[PAAt1p[i,j,k,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(1\), \(At\), \(+\)]]\)"];
DefTensor[PAate1p[i,j,k,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(Aa\), SubsuperscriptBox[\(1\), \(te\), \(+\)]]\)"];
DefTensor[PAtA1p[i,j,k,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(1\), \(tA\), \(+\)]]\)"];
DefTensor[PAt1p[i,j,k,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(1\), \(t\), \(+\)]]\)"];
(* definition for spin-2^- *)
DefTensor[PAt2m[i,j,k,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(2\), \(t\), \(-\)]]\)"];
(* definition for spin-2^+ *)
DefTensor[Pas2p[i,j,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(a\), SubsuperscriptBox[\(2\), \(s\), \(+\)]]\)"];
DefTensor[PaAst2p[i,j,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(aA\), SubsuperscriptBox[\(2\), \(st\), \(+\)]]\)"];
DefTensor[PAats2p[i,j,k,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(Aa\), SubsuperscriptBox[\(2\), \(ts\), \(+\)]]\)"];
DefTensor[PAt2p[i,j,k,c,a,b],M,PrintAs->"\!\(\*SubsuperscriptBox[\(P\), \(A\), SubsuperscriptBox[\(2\), \(t\), \(+\)]]\)"];
(* rules for spin-0^- *)
rulespin0m={PAA0m[i_,j_,k_,c_,a_,b_]:>Antisymmetrize[Antisymmetrize[2/3 \[Theta][i,c]\[Theta][j,a]\[Theta][k,b]+1/3 \[Theta][i,a]\[Theta][j,b]\[Theta][k,c],{i,j}],{a,b}]};
(* rules for spin-0^+ *)
rulespin0p={Pat0p[i_,j_,a_,b_]:>\[Omega][i,j]\[Omega][a,b],
Pats0p[i_,j_,a_,b_]:>1/Sqrt[3] \[Omega][i,j]\[Theta][a,b],
PaAtV0p[i_,j_,c_,a_,b_]:>Sqrt[2]/Sqrt[3] Antisymmetrize[k1[b]\[Omega][i,j]\[Theta][c,a],{a,b}],
Past0p[i_,j_,a_,b_]:>1/Sqrt[3] \[Theta][i,j] \[Omega][a,b],
Pas0p[i_,j_,a_,b_]:>1/3 \[Theta][i,j]\[Theta][a,b],
PaAsV0p[i_,j_,c_,a_,b_]:>Sqrt[2]/3 Antisymmetrize[k1[b]\[Theta][i,j]\[Theta][c,a],{a,b}],
PAaVt0p[i_,j_,k_,a_,b_]:>Sqrt[2]/Sqrt[3] Antisymmetrize[k1[j] \[Theta][k,i] \[Omega][b,a],{i,j}],
PAaVs0p[i_,j_,k_,a_,b_]:>Sqrt[2]/3 Antisymmetrize[k1[j] \[Theta][i,k] \[Theta][b,a],{i,j}],
PAV0p[i_,j_,k_,c_,a_,b_]:>2/3 Antisymmetrize[Antisymmetrize[\[Theta][k,j]\[Omega][i,a] \[Theta][c,b],{i,j}],{a,b}]};
(* rules for spin-1^- *)
rulespin1m={Pas1m[i_,j_,a_,b_]:>2Symmetrize[Symmetrize[\[Theta][i,a]\[Omega][j,b],{i,j}],{a,b}],
Pasm1m[i_,j_,a_,b_]:>2Antisymmetrize[Symmetrize[\[Theta][i,a]\[Omega][j,b],{i,j}],{a,b}],
PaAsV1m[i_,j_,c_,a_,b_]:>Sqrt[2]Antisymmetrize[Symmetrize[k1[j]\[Theta][i,a]\[Theta][c,b],{i,j}],{a,b}],
PaAst1m[i_,j_,c_,a_,b_]:>2Antisymmetrize[Symmetrize[k1[b]\[Theta][i,a]\[Omega][c,j],{i,j}],{a,b}],
Pams1m[i_,j_,a_,b_]:>2Symmetrize[Antisymmetrize[\[Theta][i,a]\[Omega][j,b],{i,j}],{a,b}],
Pam1m[i_,j_,a_,b_]:>2Antisymmetrize[Antisymmetrize[\[Theta][i,a]\[Omega][j,b],{i,j}],{a,b}],
PaAmV1m[i_,j_,c_,a_,b_]:>Sqrt[2]Antisymmetrize[Antisymmetrize[k1[j]\[Theta][i,a]\[Theta][c,b],{i,j}],{a,b}],
PaAmt1m[i_,j_,c_,a_,b_]:>2Antisymmetrize[Antisymmetrize[k1[b]\[Theta][i,a]\[Omega][c,j],{i,j}],{a,b}],
PAaVs1m[i_,j_,k_,a_,b_]:>Sqrt[2]Antisymmetrize[Symmetrize[k1[b]\[Theta][i,a]\[Theta][j,k],{a,b}],{i,j}],
PAaVm1m[i_,j_,k_,a_,b_]:>Sqrt[2]Antisymmetrize[Antisymmetrize[k1[b]\[Theta][i,a]\[Theta][j,k],{a,b}],{i,j}],
PAV1m[i_,j_,k_,c_,a_,b_]:>Antisymmetrize[Antisymmetrize[\[Theta][c,b]\[Theta][i,a]\[Theta][j,k],{i,j}],{a,b}],
PAVt1m[i_,j_,k_,c_,a_,b_]:>Sqrt[2]Antisymmetrize[Antisymmetrize[\[Omega][c,b]\[Theta][i,a]\[Theta][j,k],{a,b}],{i,j}],
PAats1m[i_,j_,k_,a_,b_]:>2Antisymmetrize[Symmetrize[k1[j]\[Theta][i,a]\[Omega][k,b],{a,b}],{i,j}],
PAatm1m[i_,j_,k_,a_,b_]:>2Antisymmetrize[Antisymmetrize[k1[j]\[Theta][i,a]\[Omega][k,b],{a,b}],{i,j}],
PAtV1m[i_,j_,k_,c_,a_,b_]:>Sqrt[2]Antisymmetrize[Antisymmetrize[\[Theta][c,b]\[Theta][i,a]\[Omega][j,k],{a,b}],{i,j}],
PAt1m[i_,j_,k_,c_,a_,b_]:>2Antisymmetrize[Antisymmetrize[\[Omega][c,b]\[Theta][i,a]\[Omega][j,k],{i,j}],{a,b}]};
(*rules for spin-1^+*)
rulespin1p={Pae1p[i_,j_,a_,b_]:>Antisymmetrize[Antisymmetrize[\[Theta][i,a]\[Theta][j,b],{i,j}],{a,b}],
PaAeA1p[i_,j_,c_,a_,b_]:>Sqrt[2]Antisymmetrize[Antisymmetrize[k1[b]\[Theta][i,a]\[Theta][j,c],{i,j}],{a,b}],
PaAet1p[i_,j_,c_,a_,b_]:>Antisymmetrize[Antisymmetrize[k1[c]\[Theta][i,a]\[Theta][j,b],{i,j}],{a,b}],
PAaAe1p[i_,j_,k_,a_,b_]:>Sqrt[2]Antisymmetrize[Antisymmetrize[k1[j]\[Theta][i,a]\[Theta][k,b],{a,b}],{i,j}],
PAA1p[i_,j_,k_,c_,a_,b_]:>Antisymmetrize[Antisymmetrize[\[Theta][i,c]\[Theta][k,b]\[Omega][j,a]+\[Theta][i,a]\[Theta][k,c]\[Omega][j,b],{i,j}],{a,b}],
PAAt1p[i_,j_,k_,c_,a_,b_]:>-Sqrt[2]Antisymmetrize[Antisymmetrize[\[Theta][j,a]\[Theta][k,b]\[Omega][i,c],{a,b}],{i,j}],
PAate1p[i_,j_,k_,a_,b_]:>Antisymmetrize[Antisymmetrize[k1[k]\[Theta][i,a]\[Theta][j,b],{a,b}],{i,j}],
PAtA1p[i_,j_,k_,c_,a_,b_]:>-Sqrt[2]Antisymmetrize[Antisymmetrize[\[Theta][i,b]\[Theta][j,c]\[Omega][k,a],{a,b}],{i,j}],
PAt1p[i_,j_,k_,c_,a_,b_]:>Antisymmetrize[Antisymmetrize[\[Theta][i,a]\[Theta][j,b]\[Omega][k,c],{i,j}],{a,b}]};
(*rules for spin-2^-*)
rulespin2m={PAt2m[i_,j_,k_,c_,a_,b_]:>Antisymmetrize[Antisymmetrize[2/3 \[Theta][i,c]\[Theta][j,b]\[Theta][k,a]+2/3 \[Theta][i,a]\[Theta][j,b]\[Theta][k,c]-\[Theta][b,c]\[Theta][j,k]\[Theta][i,a],{i,j}],{a,b}]};
(*rules for spin-2^+*)
rulespin2p={Pas2p[i_,j_,a_,b_]:>Symmetrize[Symmetrize[\[Theta][i,a]\[Theta][j,b]-1/3 \[Theta][i,j]\[Theta][a,b],{i,j}],{a,b}],
PaAst2p[i_,j_,c_,a_,b_]:>Sqrt[2]Antisymmetrize[Symmetrize[k1[b](\[Theta][i,a]\[Theta][j,c]-1/3 \[Theta][i,j]\[Theta][a,c]),{i,j}],{a,b}],
PAats2p[i_,j_,k_,a_,b_]:>Sqrt[2]Antisymmetrize[Symmetrize[k1[j](\[Theta][i,a]\[Theta][k,b]-1/3 \[Theta][i,k]\[Theta][a,b]),{a,b}],{i,j}],
PAt2p[i_,j_,k_,c_,a_,b_]:>Antisymmetrize[Antisymmetrize[\[Theta][i,c]\[Theta][k,a]\[Omega][j,b]+\[Theta][i,a]\[Theta][k,c]\[Omega][j,b]-2/3 \[Theta][b,c]\[Theta][k,j]\[Omega][i,a],{i,j}],{a,b}]};
(* \:603b\:89c4\:5219 *)
rulespin=Flatten@{rulespin0m,rulespin0p,rulespin1m,rulespin1p,rulespin2m,rulespin2p};
P0m[i_,j_,k_,c_,a_,b_]:=DiagonalMatrix[{0,0,0,0,PAA0m[i,j,k,c,a,b],0}];
P0p[i_,j_,k_,c_,a_,b_]:={{Pat0p[i,j,a,b],Pats0p[i,j,a,b],0,PaAtV0p[i,j,c,a,b],0,0},
{Past0p[i,j,a,b],Pas0p[i,j,a,b],0,PaAsV0p[i,j,c,a,b],0,0},
{0,0,0,0,0,0},
{PAaVt0p[i,j,k,a,b],PAaVs0p[i,j,k,a,b],0,PAV0p[i,j,k,c,a,b],0,0},
{0,0,0,0,0,0},
{0,0,0,0,0,0}};
P1m[i_,j_,k_,c_,a_,b_]:={{0,0,0,0,0,0},
{0,Pas1m[i,j,a,b],Pasm1m[i,j,a,b],PaAsV1m[i,j,c,a,b],0,PaAst1m[i,j,c,a,b]},
{0,Pams1m[i,j,a,b],Pam1m[i,j,a,b],PaAmV1m[i,j,c,a,b],0,PaAmt1m[i,j,c,a,b]},
{0,PAaVs1m[i,j,k,a,b],PAaVm1m[i,j,k,a,b],PAV1m[i,j,k,c,a,b],0,PAVt1m[i,j,k,c,a,b]},
{0,0,0,0,0,0},
{0,PAats1m[i,j,k,a,b],PAatm1m[i,j,k,a,b],PAtV1m[i,j,k,c,a,b],0,PAt1m[i,j,k,c,a,b]}};
P1p[i_,j_,k_,c_,a_,b_]:={{0,0,0,0,0,0},
{0,0,0,0,0,0},
{0,0,Pae1p[i,j,a,b],0,PaAeA1p[i,j,c,a,b],PaAet1p[i,j,c,a,b]},
{0,0,0,0,0,0},
{0,0,PAaAe1p[i,j,k,a,b],0,PAA1p[i,j,k,c,a,b],PAAt1p[i,j,k,c,a,b]},
{0,0,PAate1p[i,j,k,a,b],0,PAtA1p[i,j,k,c,a,b],PAt1p[i,j,k,c,a,b]}};
P2m[i_,j_,k_,c_,a_,b_]:=DiagonalMatrix[{0,0,0,0,0,PAt2m[i,j,k,c,a,b]}];
P2p[i_,j_,k_,c_,a_,b_]:={{0,0,0,0,0,0},
{0,Pas2p[i,j,a,b],0,0,0,PaAst2p[i,j,c,a,b]},
{0,0,0,0,0,0},
{0,0,0,0,0,0},
{0,0,0,0,0,0},
{0,PAats2p[i,j,k,a,b],0,0,0,PAt2p[i,j,k,c,a,b]}};
spinProj=(P0m[i,j,k,c,a,b]+P0p[i,j,k,c,a,b]+P1m[i,j,k,c,a,b]+P1p[i,j,k,c,a,b]+P2m[i,j,k,c,a,b]+P2p[i,j,k,c,a,b]//.rulespin);

ruleikv={kv[\[Mu]_]:>I kv[\[Mu]],kv2[]:>-kv2[]};

Oa[a_,\[Mu]_,\[Nu]_,c_]:=Module[{\[Sigma],\[Tau]},8A1 metric[ChangeIndex[\[Sigma]],ChangeIndex[\[Tau]]]Antisymmetrize[metric[a,c]metric[\[Nu],\[Mu]]kv[\[Sigma]]kv[\[Tau]],{\[Mu],\[Sigma]}]+4A2 (Antisymmetrize[metric[\[Nu],\[Mu]]kv[a]kv[c],{\[Mu],a}]+metric[ChangeIndex[\[Sigma]],ChangeIndex[\[Tau]]]Antisymmetrize[metric[c,\[Mu]]metric[a,\[Nu]]kv[\[Sigma]]kv[\[Tau]],{\[Nu],\[Sigma]}])+4A3 (metric[ChangeIndex[\[Sigma]],ChangeIndex[\[Tau]]]Antisymmetrize[metric[a,\[Mu]]metric[\[Nu],c]kv[\[Sigma]]kv[\[Tau]],{\[Mu],\[Sigma]}]-metric[ChangeIndex[\[Sigma]],ChangeIndex[\[Tau]]]Antisymmetrize[metric[a,\[Mu]]metric[\[Nu],\[Sigma]]kv[c]kv[\[Tau]],{\[Mu],\[Sigma]}])]/.ruleikv//ContractMetric;
OaA[a_,\[Mu]_,\[Nu]_,c_,d_]:=Module[{\[Sigma],\[Tau]},2b0( Antisymmetrize[metric[\[Mu],c]metric[\[Nu],d]kv[a],{c,d}]-Antisymmetrize[metric[a,\[Nu]]metric[\[Mu],c]kv[d],{c,d}]-Antisymmetrize[metric[a,\[Mu]]metric[\[Nu],d]kv[c],{d,c}])+4A1 (Antisymmetrize[metric[a,d]metric[c,\[Mu]]kv[\[Nu]],{d,c}]-Antisymmetrize[metric[a,d]metric[\[Nu],\[Mu]]kv[c],{d,c}])+2A2 (2Antisymmetrize[metric[\[Nu],a]metric[\[Mu],c]kv[d],{c,d}]-Antisymmetrize[metric[\[Nu],\[Mu]]metric[a,c]kv[d],{c,d}]+Antisymmetrize[metric[d,\[Mu]]metric[a,c]kv[\[Nu]],{c,d}])-2A3 (Antisymmetrize[metric[a,\[Mu]]metric[\[Nu],d]kv[c],{d,c}]-Antisymmetrize[metric[c,\[Mu]]metric[\[Nu],d]kv[a],{c,d}])]/.ruleikv//ContractMetric;
OAa[a_,b_,\[Mu]_,\[Nu]_,c_]:=Module[{\[Sigma],\[Tau]},-2b0(Antisymmetrize[metric[\[Mu],c]metric[\[Nu],b]kv[a],{b,a}]+Antisymmetrize[metric[\[Mu],a]metric[\[Nu],c]kv[b],{a,b}]-Antisymmetrize[metric[\[Mu],a]metric[\[Nu],b]kv[c],{a,b}])-4A1 (Antisymmetrize[metric[\[Nu],a]metric[b,c]kv[\[Mu]],{a,\[Mu]}]-Antisymmetrize[metric[\[Nu],b]metric[a,c]kv[\[Mu]],{b,\[Mu]}])-2A2 (2Antisymmetrize[metric[\[Mu],c]metric[\[Nu],a]kv[b],{a,b}]+Antisymmetrize[metric[b,c]metric[\[Nu],\[Mu]]kv[a],{\[Mu],a}]-Antisymmetrize[metric[a,c]metric[\[Nu],\[Mu]]kv[b],{\[Mu],b}])-2A3 (Antisymmetrize[metric[\[Nu],c]metric[\[Mu],a]kv[b],{c,b}]-Antisymmetrize[metric[\[Nu],c]metric[\[Mu],b]kv[a],{c,a}])]/.ruleikv//ContractMetric;
OA[a_,b_,\[Mu]_,\[Nu]_,c_,d_]:=Module[{\[Sigma],\[Tau]},-2b0 (Antisymmetrize[Antisymmetrize[metric[\[Mu],d]metric[\[Nu],a]metric[c,b],{a,b}],{c,d}]-Antisymmetrize[Antisymmetrize[metric[\[Mu],a]metric[\[Nu],d]metric[c,b],{a,b}],{c,d}])-4A1 (Antisymmetrize[Antisymmetrize[metric[\[Nu],\[Mu]]metric[a,c]metric[b,d],{\[Mu],a}],{c,d}]-Antisymmetrize[Antisymmetrize[metric[\[Nu],\[Mu]]metric[b,c]metric[a,d],{\[Mu],b}],{c,d}])-2A2 (2Antisymmetrize[Antisymmetrize[metric[\[Nu],b]metric[a,c]metric[\[Mu],d],{b,a}],{c,d}]+Antisymmetrize[Antisymmetrize[metric[\[Nu],a]metric[\[Mu],c]metric[b,d],{a,\[Mu]}],{c,d}]-Antisymmetrize[Antisymmetrize[metric[\[Nu],b]metric[\[Mu],c]metric[a,d],{b,\[Mu]}],{c,d}])+2A3 Antisymmetrize[Antisymmetrize[metric[\[Nu],d]metric[\[Mu],a]metric[c,b],{a,b}],{c,d}]+8B1 metric[ChangeIndex[\[Sigma]],ChangeIndex[\[Tau]]]Antisymmetrize[Antisymmetrize[metric[a,c]metric[b,d]metric[\[Nu],\[Mu]]kv[\[Sigma]]kv[\[Tau]],{\[Mu],\[Sigma]}],{c,d}]+8B2 Antisymmetrize[Antisymmetrize[metric[\[Mu],d]metric[\[Nu],b]kv[c]kv[a],{b,a}],{d,c}]+2B3 (Antisymmetrize[Antisymmetrize[metric[\[Mu],a]metric[c,b]kv[d]kv[\[Nu]],{a,b}],{c,d}]-metric[ChangeIndex[\[Sigma]],ChangeIndex[\[Tau]]]Antisymmetrize[Antisymmetrize[metric[\[Mu],a]metric[c,b]metric[d,\[Nu]]kv[\[Sigma]]kv[\[Tau]],{a,b}],{c,d}]+Antisymmetrize[Antisymmetrize[metric[d,a]metric[c,\[Nu]]kv[b]kv[\[Mu]],{a,b}],{c,d}]+Antisymmetrize[Antisymmetrize[metric[c,a]metric[\[Mu],\[Nu]]kv[b]kv[d],{a,b}],{c,d}])+2B4 (Antisymmetrize[Antisymmetrize[metric[\[Mu],a]metric[c,\[Nu]]kv[b]kv[d],{a,b}],{c,d}]+Antisymmetrize[Antisymmetrize[metric[c,\[Mu]]metric[a,\[Nu]]kv[b]kv[d],{a,b}],{c,d}])+2B5 (metric[ChangeIndex[\[Sigma]],ChangeIndex[\[Tau]]]Antisymmetrize[Antisymmetrize[metric[\[Mu],c]metric[\[Nu],a]metric[b,d]kv[\[Sigma]]kv[\[Tau]],{a,b}],{d,c}]+Antisymmetrize[Antisymmetrize[metric[\[Mu],c]metric[a,d]kv[b]kv[\[Nu]],{a,b}],{c,d}]+Antisymmetrize[metric[\[Nu],b]metric[a,d]kv[\[Mu]]kv[c],{b,\[Mu]}]-Antisymmetrize[metric[\[Nu],a]metric[b,d]kv[\[Mu]]kv[c],{a,\[Mu]}])]/.ruleikv//ContractMetric;
Oo[a_,b_,\[Mu]_,\[Nu]_,c_,d_]:={{Symmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}],Symmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}],Antisymmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}],Antisymmetrize[Symmetrize[OaA[a,b,\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Symmetrize[OaA[a,b,\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Symmetrize[OaA[a,b,\[Nu],c,d],{a,b}],{c,d}]},{Symmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}],Symmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}],Antisymmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}],Antisymmetrize[Symmetrize[OaA[a,b,\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Symmetrize[OaA[a,b,\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Symmetrize[OaA[a,b,\[Nu],c,d],{a,b}],{c,d}]},{Symmetrize[Antisymmetrize[Oa[a,b,c,d],{a,b}],{c,d}],Symmetrize[Antisymmetrize[Oa[a,b,c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[Oa[a,b,c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OaA[a,b,\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OaA[a,b,\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OaA[a,b,\[Nu],c,d],{a,b}],{c,d}]},{Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}],Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Nu],c,d],{a,b}],{c,d}]},{Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}],Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Nu],c,d],{a,b}],{c,d}]},{Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}],Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Nu],c,d],{a,b}],{c,d}],Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Nu],c,d],{a,b}],{c,d}]}}
S[a_,b_,c_,\[Gamma]_,\[Alpha]_,\[Beta]_]:={{1,0,0},{0,1,0},{Symmetrize[kv[a]metric[\[Alpha],c]metric[\[Beta],b]-kv[b]metric[\[Alpha],c]metric[\[Beta],a],{\[Alpha],\[Beta]}],-kv[c]Antisymmetrize[metric[\[Alpha],a]metric[\[Beta],b],{\[Alpha],\[Beta]}],1/2 Antisymmetrize[metric[\[Alpha],c] metric[\[Beta],a]metric[\[Gamma],b]+metric[\[Alpha],b] metric[\[Beta],a]metric[\[Gamma],c]+metric[\[Alpha],b] metric[\[Beta],c]metric[\[Gamma],a],{\[Alpha],\[Beta]}]}}/.ruleikv;
Os[a_,b_,\[Mu]_,\[Nu]_,c_,d_]:={{Symmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Symmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Antisymmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[2]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]]},{Symmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Symmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Antisymmetrize[Symmetrize[Oa[a,b,c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[2]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Symmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}]S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]]},{Symmetrize[Antisymmetrize[Oa[a,b,c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Symmetrize[Antisymmetrize[Oa[a,b,c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Antisymmetrize[Antisymmetrize[Oa[a,b,c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[2]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OaA[a,b,\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]]},{Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Antisymmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[2]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]]},{Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Antisymmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[2]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]]},{Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Symmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[1]]],Antisymmetrize[Antisymmetrize[OAa[a,b,\[Mu],c,d],{a,b}],{c,d}]+Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[2]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]],Module[{\[Alpha],\[Beta],\[Gamma]},Antisymmetrize[Antisymmetrize[OA[a,b,\[Mu],\[Gamma],\[Alpha],\[Beta]],{a,b}],{\[Alpha],\[Beta]}] S[ChangeIndex[\[Alpha]],ChangeIndex[\[Beta]],ChangeIndex[\[Gamma]],\[Nu],c,d][[3]][[3]]]}}
(* GR *)
Ogr[\[Mu]_,\[Nu]_,\[Rho]_,\[Sigma]_]:=(Symmetrize[metric[\[Mu],\[Rho]]metric[\[Sigma],\[Nu]],{\[Rho],\[Sigma]}]-metric[\[Mu],\[Nu]]metric[\[Rho],\[Sigma]])kv2[]+metric[\[Mu],\[Nu]]kv[\[Rho]]kv[\[Sigma]]+metric[\[Rho],\[Sigma]]kv[\[Mu]]kv[\[Nu]]-(Symmetrize[metric[\[Nu],\[Rho]]kv[\[Sigma]]kv[\[Mu]]+metric[\[Mu],\[Rho]]kv[\[Sigma]]kv[\[Nu]],{\[Rho],\[Sigma]}]);
OGR[\[Mu]_,\[Nu]_,\[Rho]_,\[Sigma]_]:=Ogr[\[Mu],\[Nu],\[Rho],\[Sigma]]{{1,1},{1,1}};
Pgr0p[i_,j_,a_,b_]:={{Pat0p[i,j,a,b],Pats0p[i,j,a,b]},
{Past0p[i,j,a,b],Pas0p[i,j,a,b]}};
Pgr1m[i_,j_,a_,b_]:={{0,0},
{0,Pas1m[i,j,a,b]}};
Pgr2p[i_,j_,a_,b_]:={{0,0},
{0,Pas2p[i,j,a,b]}};
(* parameters *)
ruleparas={A1->a1/2-a3/18,A2->a1/2+a3/9,A3->-a1/2+a2/9,B1->b1/6+b3/2+b5/4-b6,B2->b1/6+b2/4-b3/2,B3->b2/4-b3+b4/2-b5/2+4b6,B4->-3 b2/4+b3-b4/2+b5/2,B5->-2 b1/3+b2/2};
(*ruleparas1={b0->\[Lambda],A1->1/12 (4a+b+3\[Lambda]),A2->-(1/6)(-2a+b-3\[Lambda]),A3->1/3 (-a+2c-3\[Lambda]),B1->1/6 (2p+q),B2->1/6 (2p+q-6r),B5->2/3 (p-q),B3->s+t,B4->s-t};*)
rulegr={A1->0,A2->0,A3->0,B1->0,B2->0,B3->0,B4->0,B5->0,b0->1};
(*invruleparam1=Solve[Table[ruleparas1[[ii]][[1]]==ruleparas1[[ii]][[2]],{ii,1,Length@ruleparas1}],{\[Lambda],a,b,c,p,q,r,s,t}][[1]];*)


PPo={P0m[i,j,k,c,a,b],P0p[i,j,k,c,a,b],P1m[i,j,k,c,a,b],P1p[i,j,k,c,a,b],P2m[i,j,k,c,a,b],P2p[i,j,k,c,a,b]};
PPt={P0m[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]],P0p[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]],P1m[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]],P1p[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]],P2m[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]],P2p[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]]};
PPp={P0m[i,j,k,\[Tau],\[Xi],\[Zeta]],P0p[i,j,k,\[Tau],\[Xi],\[Zeta]],P1m[i,j,k,\[Tau],\[Xi],\[Zeta]],P1p[i,j,k,\[Tau],\[Xi],\[Zeta]],P2m[i,j,k,\[Tau],\[Xi],\[Zeta]],P2p[i,j,k,\[Tau],\[Xi],\[Zeta]]};
POP=Table[Table[Module[{xxx,yyy},xxx=PPo[[iii]][[ii]][[ii]]Oo[-a,-b,-c,-\[Sigma],-\[Mu],-\[Nu]][[ii]][[jj]]PPt[[iii]][[jj]][[jj]]//.rulespin//.ruleProj//ToCanonical//ContractMetric//Simplify;yyy=PPp[[iii]][[ii]][[jj]]//.rulespin//.ruleProj//ToCanonical//Simplify;If[NumberQ@xxx==True,0,Collect[-1/2 xxx/yyy//Simplify,kv2[]]]],{ii,1,6},{jj,1,6}],{iii,1,6}];


PPo={P0m[i,j,k,c,a,b],P0p[i,j,k,c,a,b],P1m[i,j,k,c,a,b],P1p[i,j,k,c,a,b],P2m[i,j,k,c,a,b],P2p[i,j,k,c,a,b]};
PPt={P0m[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]],P0p[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]],P1m[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]],P1p[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]],P2m[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]],P2p[\[Mu],\[Nu],\[Sigma],\[Tau],\[Xi],\[Zeta]]};
PPp={P0m[i,j,k,\[Tau],\[Xi],\[Zeta]],P0p[i,j,k,\[Tau],\[Xi],\[Zeta]],P1m[i,j,k,\[Tau],\[Xi],\[Zeta]],P1p[i,j,k,\[Tau],\[Xi],\[Zeta]],P2m[i,j,k,\[Tau],\[Xi],\[Zeta]],P2p[i,j,k,\[Tau],\[Xi],\[Zeta]]};
POsP=Table[Table[Module[{xxx,yyy},xxx=PPo[[iii]][[ii]][[ii]]Os[-a,-b,-c,-\[Sigma],-\[Mu],-\[Nu]][[ii]][[jj]]PPt[[iii]][[jj]][[jj]]//.rulespin//.ruleProj//ToCanonical//ContractMetric//Simplify;yyy=PPp[[iii]][[ii]][[jj]]//.rulespin//.ruleProj//ToCanonical//Simplify;If[NumberQ@xxx==True,0,Collect[-1/2 xxx/yyy//Simplify,kv2[]]]],{ii,1,6},{jj,1,6}],{iii,1,6}];


(*End[]*)


EndPackage[]
