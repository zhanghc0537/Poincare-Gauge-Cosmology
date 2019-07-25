(* ::Package:: *)

(************************ 0. Info and copyright ***********************)


PGC121`$Version={"1.2.1",{2019,7,22}}


(* Title: PGC, Poincare Gauge Cosmology *)

(* Author: Hongchao Zhang *)

(* Summary: Symbolic computing package on Poincare Gauge Cosmology *)

(* Copyright (C) 2018-2019 Hongchao Zhang *)

(* Email: zhanghc@mail.dlut.edu.cn *)

(* Affiliation: Dalian University of Technology & Penn State University *)

(* Package Version: 1.2.1 *)

(* Brief Discussion: 
   - Define a 4D metric-affinely connected manifold.
   - Define a Chart and a time-like observer.
   - Define the components of metric and torsion tensor on the FLRW cosmology
   - Vary the given gravitational Lagrangian with respect the metric and torsion, to get the Einstein and the Cartan field equation, respectively
   - Calculate the cosmological equations on the FLRW background.
   - Export and store the results in files in the working directory.
*)

(* Users are required to cite the PGC papers: *)

(* arXiv: 1904.03545, arXiv: 1906.04340 *)

(* This package is based on the xAct series, version 1.1.2. For more information please transfer to the site: http://www.xact.es/ *)


(************************ 1. Begin package ***********************)


BeginPackage["PGC121`",{"xAct`xCore`","xAct`xPerm`","xAct`xTensor`","xAct`xCoba`","xAct`xPert`"}]


(* Parallelization setting up *)
ParallelNeeds["xAct`xCore`"];
ParallelNeeds["xAct`xPerm`"];
ParallelNeeds["xAct`xTensor`"];
ParallelNeeds["xAct`xCoba`"];
ParallelNeeds["xAct`xPert`"];


PGC121`$barslength=60;
PGC121`bars=StringJoin[Table["-",{PGC121`$barslength}]];


Print[PGC121`bars];
Print[PGC121`bars];
Print["Package PGC`  version ",PGC121`$Version[[1]],", ",PGC121`$Version[[2]]];
(*Print["PGC, Poincare Gauge Cosmology"];
Print["Symbolic computing package on Poincare Gauge Cosmology"];
Print["Copyright (C) 2018-2019 Hongchao Zhang"];
Print["Email: zhanghc@mail.dlut.edu.cn"];
Print["Dalian University of Technology & Penn State University"];
Print["This program is based on the xAct series, version 1.1.2: http://www.xact.es/"];
Print["Users are required to cite the PGC papers:"];
Print["arXiv: 1904.03545, arXiv: 1906.04340"];
Print["Version 1.2.1: Nine-parameter gravitational Lagrangian, Einstein and Cartan equations, and their components on background"];*)
Print["Title: PGC, Poincare Gauge Cosmology"];
Print["Author: Hongchao Zhang"];
Print["Summary: Symbolic computing package on Poincare Gauge Cosmology"];
Print["Copyright (C) 2018-2019 Hongchao Zhang"];
Print["Email: zhanghc@mail.dlut.edu.cn"];
Print["Affiliation: Dalian University of Technology & Penn State University"];
Print["Package Version: 1.2.1"];
Print["Brief Discussion: 
   - Define a 4D metric-affinely connected manifold.
   - Define a Chart and a time-like observer.
   - Define the components of metric and torsion tensor on the FLRW cosmology
   - Vary the given gravitational Lagrangian with respect the metric and torsion, to get the Einstein and the Cartan field equation, respectively
   - Calculate the cosmological equations on the FLRW background.
   - Export and store the results in files in the working directory."];
Print["Users are required to cite the PGC papers:"];
Print["arXiv: 1904.03545, arXiv: 1906.04340"];
Print["This package is based on the xAct series, version 1.1.2. For more information please transfer to the site: http://www.xact.es/"];
Print[PGC121`bars];


(*usage*)
EinsteinFieldEq::usage="EinsteinFieldEq[var_List] exports the gravitational part of the Einstein field equation corresponding to the Lagrangian assigned by a List of parameters.";
CartanFieldEq::usage="CartanFieldEq[var_List] exports the gravitational part of the Cartan field equation corresponding to the Lagrangian assigned by a List of parameters.";
EinsteinCompEqDim::usage="EinsteinCompEqDim[\[Mu]_Integer,\[Nu]_Integer,var_List] exports the {\[Mu],\[Nu]} (\[Mu],\[Nu]=0,1,2,3) component of the Einstein field equation in the FLRW case, corresponding to the Lagrangian assigned by a List of parameters.";
CartanCompEqDim::usage="CartanCompEqDim[\[Mu]_Integer,\[Nu]_Integer,\[Rho]_Integer,var_List] exports the {\[Mu],\[Nu],\[Rho]} (\[Mu],\[Nu],\[Rho]=0,1,2,3) component of the Cartan field equation in the FLRW case, corresponding to the Lagrangian assigned by a List of parameters.";
EMTConservedCompEqDim::usage="EMTConservedCompEqDim[\[Mu]_Integer] exports the {\[Mu]} (\[Mu]=0,1,2,3) component of the conservation equation of energy-momentum tensor in the FLRW case.";
LagrangianCompvar::usage="LagrangianCompvar[var_List] exports the component of the Lagrangian assigned by a List of parameters in the FLRW case.";


(*Begin["`Private`"]*)


(************************ 2. General Statements ***********************)


$DefInfoQ=False;
$CVVerbose=False;
$PrePrint=ScreenDollarIndices;


(* Adding the working directory *)
$PackageDirectory=ParentDirectory[DirectoryName[FindFile["PGC121`"]]];
If[!ContainsAny[$Path,{$PackageDirectory}],PrependTo[$Path,$PackageDirectory]];
Print["The current working directory is: ",$PackageDirectory];
Print[PGC121`bars];
PrintTemporary["Wait a moment, the code is running ..."];
Print[PGC121`bars];


(* my options: Torsion or Contorsion *)
$WhichVar=1;
(* Pre evaluating cdcdTorsionCD or not: 1 is yes *)
$PreEvaluatecdcdTorsionCD=1;


(* Define a manifold M *)
DefManifold[M,4,{\[Mu],\[Nu],\[Sigma],\[Tau],\[Theta],\[Upsilon],\[Xi],\[Zeta],\[Beta],\[Gamma]},PrintAs->"\!\(\*SuperscriptBox[\(M\), \(4\)]\)"];


(* Define a metric g and its perturbation *)
DefMetric[-1,metric[-\[Mu],-\[Nu]],cd,SymbolOfCovD->{";","D"},WeightedWithBasis->AIndex,PrintAs->"g"];
DefMetricPerturbation[metric,metpert,\[Epsilon]];
PrintAs[metpert]^="\[Delta]g";


Unprotect[IndexForm];
IndexForm[LI[x_]]:=ColorString[ToString[x],RGBColor[0,0,1]];
Protect[IndexForm];


(* Define a metric-affine connection compatible with the general covariant derivative CD and metric g *)
DefCovD[CD[-\[Mu]],Torsion->True,SymbolOfCovD->{"|","\[Del]"},FromMetric->metric];


(* Rules among \[CapitalGamma], T and K *)
rule\[CapitalGamma]=ChristoffelCD[\[Sigma]_,\[Mu]_,\[Nu]_]->Christoffelcd[\[Sigma],\[Mu],\[Nu]]+Contorsion[\[Sigma],\[Mu],\[Nu]];
ruleK=Contorsion[\[Sigma]_,\[Mu]_,\[Nu]_]->1/2 (TorsionCD[\[Sigma],\[Mu],\[Nu]]+TorsionCD[\[Mu],\[Sigma],\[Nu]]+TorsionCD[\[Nu],\[Sigma],\[Mu]]);
(*rule=Which[$WhichVar==1,rule\[CapitalGamma]/.ruleK,$WhichVar==2,rule\[CapitalGamma]];*)


metorder=metpert[LI[order_],__]:>0/;order>1;


Clear[\[Alpha],a1,a2,a3,b1,b2,b3,b4,b5,b6];
Clear[c1];
DefConstantSymbol[\[Alpha]];
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
DefConstantSymbol[kappa,PrintAs->"\[Kappa]"];


(* Define a chart *)
DefChart[chart,M,{0,1,2,3},{t[],x[],y[],z[]},ChartColor->RGBColor[1,0,0]]


(* the background values of metric on this chart, where \[Tau] denotes the cosmic time *)
metricbg={
 {-1, 0, 0, 0},
 {0, a[t[]]^2, 0, 0},
 {0, 0, a[t[]]^2, 0},
 {0, 0, 0, a[t[]]^2}
};
(* Write the values into the background metric *)
AllComponentValues[Perturbation[metric[-\[Mu],-\[Nu]],0]//ToBasis[chart],metricbg];
MetricCompute[metric,chart,All,CVSimplify->Simplify];


metricP1=a[t[]]^2 {{-(2/a[t[]]^2) \[Phi][t[],x[],y[],z[]],2/a[t[]] PDchart[{1,-chart}][B[t[],x[],y[],z[]]],2/a[t[]] PDchart[{2,-chart}][B[t[],x[],y[],z[]]],2/a[t[]] PDchart[{3,-chart}][B[t[],x[],y[],z[]]]},{2/a[t[]] PDchart[{1,-chart}][B[t[],x[],y[],z[]]],-2 \[Psi][t[],x[],y[],z[]]+2 PDchart[{1,-chart}][PDchart[{1,-chart}][\[CapitalEpsilon][t[],x[],y[],z[]]]],2 PDchart[{1,-chart}][PDchart[{2,-chart}][\[CapitalEpsilon][t[],x[],y[],z[]]]],2 PDchart[{1,-chart}][PDchart[{3,-chart}][\[CapitalEpsilon][t[],x[],y[],z[]]]]},{2/a[t[]] PDchart[{2,-chart}][B[t[],x[],y[],z[]]],2 PDchart[{2,-chart}][PDchart[{1,-chart}][\[CapitalEpsilon][t[],x[],y[],z[]]]],-2 \[Psi][t[],x[],y[],z[]]+2 PDchart[{2,-chart}][PDchart[{2,-chart}][\[CapitalEpsilon][t[],x[],y[],z[]]]],2 PDchart[{2,-chart}][PDchart[{3,-chart}][\[CapitalEpsilon][t[],x[],y[],z[]]]]},{2/a[t[]] PDchart[{3,-chart}][B[t[],x[],y[],z[]]],2 PDchart[{3,-chart}][PDchart[{1,-chart}][\[CapitalEpsilon][t[],x[],y[],z[]]]],2 PDchart[{3,-chart}][PDchart[{2,-chart}][\[CapitalEpsilon][t[],x[],y[],z[]]]],-2 \[Psi][t[],x[],y[],z[]]+2 PDchart[{3,-chart}][PDchart[{3,-chart}][\[CapitalEpsilon][t[],x[],y[],z[]]]]}};
AllComponentValues[Perturbation[metric[-\[Mu],-\[Nu]],1]//ToBasis[chart],metricP1];
metricP2=a[t[]]^2 {{-(2/a[t[]]^2) \[Phi]2[t[],x[],y[],z[]],2/a[t[]] PDchart[{1,-chart}][B2[t[],x[],y[],z[]]],2/a[t[]] PDchart[{2,-chart}][B2[t[],x[],y[],z[]]],2/a[t[]] PDchart[{3,-chart}][B2[t[],x[],y[],z[]]]},{2/a[t[]] PDchart[{1,-chart}][B2[t[],x[],y[],z[]]],-2 \[Psi]2[t[],x[],y[],z[]]+2 PDchart[{1,-chart}][PDchart[{1,-chart}][\[CapitalEpsilon]2[t[],x[],y[],z[]]]],2 PDchart[{1,-chart}][PDchart[{2,-chart}][\[CapitalEpsilon]2[t[],x[],y[],z[]]]],2 PDchart[{1,-chart}][PDchart[{3,-chart}][\[CapitalEpsilon]2[t[],x[],y[],z[]]]]},{2/a[t[]] PDchart[{2,-chart}][B2[t[],x[],y[],z[]]],2 PDchart[{2,-chart}][PDchart[{1,-chart}][\[CapitalEpsilon]2[t[],x[],y[],z[]]]],-2 \[Psi]2[t[],x[],y[],z[]]+2 PDchart[{2,-chart}][PDchart[{2,-chart}][\[CapitalEpsilon]2[t[],x[],y[],z[]]]],2 PDchart[{2,-chart}][PDchart[{3,-chart}][\[CapitalEpsilon]2[t[],x[],y[],z[]]]]},{2/a[t[]] PDchart[{3,-chart}][B2[t[],x[],y[],z[]]],2 PDchart[{3,-chart}][PDchart[{1,-chart}][\[CapitalEpsilon]2[t[],x[],y[],z[]]]],2 PDchart[{3,-chart}][PDchart[{2,-chart}][\[CapitalEpsilon]2[t[],x[],y[],z[]]]],-2 \[Psi]2[t[],x[],y[],z[]]+2 PDchart[{3,-chart}][PDchart[{3,-chart}][\[CapitalEpsilon]2[t[],x[],y[],z[]]]]}};
AllComponentValues[Perturbation[metric[-\[Mu],-\[Nu]],2]//ToBasis[chart],metricP2];


(*DefTensorPerturbation[PertChristoffelcdPDchart[LI[order],\[Xi],-\[Mu],-\[Nu]],ChristoffelcdPDchart[\[Xi],-\[Mu],-\[Nu]],M];*)


(* Torsion and its decompositions *)
torsionvaluebg={{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,-a[t[]]^2h[t[]],0,0},{a[t[]]^2 h[t[]],0,0,0},{0,0,0,2a[t[]]^3 f[t[]]},{0,0,-2 a[t[]]^3 f[t[]],0}},{{0,0,-a[t[]]^2h[t[]],0},{0,0,0,-2 a[t[]]^3 f[t[]]},{a[t[]]^2 h[t[]],0,0,0},{0,2a[t[]]^3 f[t[]],0,0}},{{0,0,0,-a[t[]]^2h[t[]]},{0,0,2a[t[]]^3 f[t[]],0},{0,-2 a[t[]]^3 f[t[]],0,0},{a[t[]]^2 h[t[]],0,0,0}}};
AllComponentValues[Perturbation[TorsionCD[-\[Sigma],-\[Mu],-\[Nu]],0]//ToBasis[chart],torsionvaluebg];
DefTensor[VTorsion[\[Sigma]],M,PrintAs->"V"];
DefTensorPerturbation[PertVTorsion[LI[order],\[Sigma]],VTorsion[\[Sigma]],M,PrintAs->"\[Delta]V"];
DefTensor[ATorsion[\[Sigma]],M,PrintAs->"A"];
DefTensorPerturbation[PertATorsion[LI[order],\[Sigma]],ATorsion[\[Sigma]],M,PrintAs->"\[Delta]A"];
VTorsionP1={\[CapitalXi][t[],x[],y[],z[]],1/a[t[]]PDchart[{1,-chart}][\[CapitalTheta][t[],x[],y[],z[]]],1/a[t[]]PDchart[{2,-chart}][\[CapitalTheta][t[],x[],y[],z[]]],1/a[t[]]PDchart[{3,-chart}][\[CapitalTheta][t[],x[],y[],z[]]]};
ATorsionP1=-2{\[CapitalPhi][t[],x[],y[],z[]],1/a[t[]]PDchart[{1,-chart}][\[CapitalPsi][t[],x[],y[],z[]]],1/a[t[]]PDchart[{2,-chart}][\[CapitalPsi][t[],x[],y[],z[]]],1/a[t[]]PDchart[{3,-chart}][\[CapitalPsi][t[],x[],y[],z[]]]};
AllComponentValues[Perturbation[VTorsion[\[Sigma]],0]//ToBasis[chart],Perturbation[1/3 TorsionCD[\[Mu],\[Sigma],-\[Mu]],0]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues];
AllComponentValues[Perturbation[ATorsion[\[Sigma]],0]//ToBasis[chart],Simplify[(Perturbation[1/6 metric[\[Sigma],\[Zeta]]epsilonmetric[-\[Zeta],-\[Mu],-\[Nu],-\[Xi]]TorsionCD[\[Mu],\[Nu],\[Xi]],0]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToCanonical)/.epsilonToetaDown[metric,chart]//ToValues,a[t[]]>0]];
AllComponentValues[Perturbation[VTorsion[\[Sigma]],1]//ToBasis[chart],VTorsionP1];
AllComponentValues[Perturbation[ATorsion[\[Sigma]],1]//ToBasis[chart],ATorsionP1];
AllComponentValues[Perturbation[TorsionCD[\[Sigma],-\[Mu],-\[Nu]],1]//ToBasis[chart],Simplify[(Perturbation[(metric[-\[Mu],-\[Xi]]VTorsion[\[Xi]]delta[\[Sigma],-\[Nu]]-metric[-\[Nu],-\[Zeta]]VTorsion[\[Zeta]]delta[\[Sigma],-\[Mu]])+metric[\[Sigma],\[Xi]]epsilonmetric[-\[Xi],-\[Mu],-\[Nu],-\[Zeta]]ATorsion[\[Zeta]],1]//ExpandPerturbation//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToCanonical)/.epsilonToetaDown[metric,chart]//ToValues,a[t[]]>0]//ToCanonical];
VTorsionP2={\[CapitalXi]2[t[],x[],y[],z[]],1/a[t[]]PDchart[{1,-chart}][\[CapitalTheta]2[t[],x[],y[],z[]]],1/a[t[]]PDchart[{2,-chart}][\[CapitalTheta]2[t[],x[],y[],z[]]],1/a[t[]]PDchart[{3,-chart}][\[CapitalTheta]2[t[],x[],y[],z[]]]};
ATorsionP2=-2{\[CapitalPhi]2[t[],x[],y[],z[]],1/a[t[]]PDchart[{1,-chart}][\[CapitalPsi]2[t[],x[],y[],z[]]],1/a[t[]]PDchart[{2,-chart}][\[CapitalPsi]2[t[],x[],y[],z[]]],1/a[t[]]PDchart[{3,-chart}][\[CapitalPsi]2[t[],x[],y[],z[]]]};
AllComponentValues[Perturbation[VTorsion[\[Sigma]],2]//ToBasis[chart],VTorsionP2];
AllComponentValues[Perturbation[ATorsion[\[Sigma]],2]//ToBasis[chart],ATorsionP2];
AllComponentValues[Perturbation[TorsionCD[\[Sigma],-\[Mu],-\[Nu]],2]//ToBasis[chart],Simplify[(Perturbation[(metric[-\[Mu],-\[Xi]]VTorsion[\[Xi]]delta[\[Sigma],-\[Nu]]-metric[-\[Nu],-\[Zeta]]VTorsion[\[Zeta]]delta[\[Sigma],-\[Mu]])+metric[\[Sigma],\[Xi]]epsilonmetric[-\[Xi],-\[Mu],-\[Nu],-\[Zeta]]ATorsion[\[Zeta]],2]//ExpandPerturbation//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToCanonical)/.epsilonToetaDown[metric,chart]//ToValues,a[t[]]>0]//ToCanonical];


(* Define an observer *)
DefTensor[u[\[Mu]],M];
u/:u[\[Mu]_]u[-\[Mu]_]:=-1;
DefTensorPerturbation[Pertu[LI[order],\[Mu]],u[\[Mu]],M];
PrintAs[Pertu]^="\[Delta]u";
AllComponentValues[Perturbation[u[\[Mu]],0]//ToBasis[chart],{1,0,0,0}];
uupP1=1/a[t[]] {-a[t[]] \[Phi][t[],x[],y[],z[]],PDchart[{1,-chart}][v[t[],x[],y[],z[]]],PDchart[{2,-chart}][v[t[],x[],y[],z[]]],PDchart[{3,-chart}][v[t[],x[],y[],z[]]]};
AllComponentValues[Perturbation[u[\[Mu]],1]//ToBasis[chart],uupP1];
delta2u0=3 \[Phi][t[],x[],y[],z[]]^2- \[Phi]2[t[],x[],y[],z[]]+4(PDchart[{1,-chart}][B[t[],x[],y[],z[]]]PDchart[{1,-chart}][v[t[],x[],y[],z[]]]+PDchart[{2,-chart}][B[t[],x[],y[],z[]]]PDchart[{2,-chart}][v[t[],x[],y[],z[]]]+PDchart[{3,-chart}][B[t[],x[],y[],z[]]]PDchart[{3,-chart}][v[t[],x[],y[],z[]]])+(PDchart[{1,-chart}][v[t[],x[],y[],z[]]]PDchart[{1,-chart}][v[t[],x[],y[],z[]]]+PDchart[{2,-chart}][v[t[],x[],y[],z[]]]PDchart[{2,-chart}][v[t[],x[],y[],z[]]]+PDchart[{3,-chart}][v[t[],x[],y[],z[]]]PDchart[{3,-chart}][v[t[],x[],y[],z[]]]);
uupP2=1/a[t[]] {a[t[]] delta2u0,PDchart[{1,-chart}][v2[t[],x[],y[],z[]]],PDchart[{2,-chart}][v2[t[],x[],y[],z[]]],PDchart[{3,-chart}][v2[t[],x[],y[],z[]]]};
AllComponentValues[Perturbation[u[\[Mu]],2]//ToBasis[chart],uupP2];


(* Define ADM (3+1) decomposition *)
Off[DefMetric::old];
DefMetric[1,metrich[-\[Mu],-\[Nu]],cdh,{":","\[DifferentialD]"},InducedFrom->{metric,u},PrintAs->"h"];
On[DefMetric::old];
DefTensorPerturbation[Pertmetrich[LI[order],-\[Mu],-\[Nu]],metrich[-\[Mu],-\[Nu]],M,PrintAs->"\[Delta]h"];
DefScalarFunction[lapse];PrintAs[lapse]^="N";
DefTensor[shift[\[Xi]],M,OrthogonalTo->{u[-\[Xi]]},PrintAs->"N"];


(* Define the perfect fluid energy-momentum tensor *)
DefTensor[EMT[\[Mu],-\[Nu]],M,PrintAs->"T"];
DefTensorPerturbation[pertEMT[LI[order],\[Mu],-\[Nu]],EMT[\[Mu],-\[Nu]],M,PrintAs->"\[Delta]T"]
DefScalarFunction[\[Delta]\[Rho]];
DefScalarFunction[\[Delta]p];
DefScalarFunction[\[CapitalPi]];
DefTensor[stress[\[Mu],-\[Nu]],M,PrintAs->"\[CapitalPi]"];
DefTensorPerturbation[pertstress[LI[order],\[Mu],-\[Nu]],stress[\[Mu],-\[Nu]],M,PrintAs->"\[Delta]\[CapitalPi]"];
stress1={{0,0,0,0},{0,PDchart[{1,-chart}][PDchart[{1,-chart}][\[CapitalPi][t[],x[],y[],z[]]]],PDchart[{1,-chart}][PDchart[{2,-chart}][\[CapitalPi][t[],x[],y[],z[]]]],PDchart[{1,-chart}][PDchart[{3,-chart}][\[CapitalPi][t[],x[],y[],z[]]]]},{0,PDchart[{2,-chart}][PDchart[{1,-chart}][\[CapitalPi][t[],x[],y[],z[]]]],PDchart[{2,-chart}][PDchart[{2,-chart}][\[CapitalPi][t[],x[],y[],z[]]]],PDchart[{2,-chart}][PDchart[{3,-chart}][\[CapitalPi][t[],x[],y[],z[]]]]},{0,PDchart[{3,-chart}][PDchart[{1,-chart}][\[CapitalPi][t[],x[],y[],z[]]]],PDchart[{3,-chart}][PDchart[{2,-chart}][\[CapitalPi][t[],x[],y[],z[]]]],PDchart[{3,-chart}][PDchart[{3,-chart}][\[CapitalPi][t[],x[],y[],z[]]]]}}-1/3 (PDchart[{1,-chart}][PDchart[{1,-chart}][\[CapitalPi][t[],x[],y[],z[]]]]+PDchart[{2,-chart}][PDchart[{2,-chart}][\[CapitalPi][t[],x[],y[],z[]]]]+PDchart[{3,-chart}][PDchart[{3,-chart}][\[CapitalPi][t[],x[],y[],z[]]]]){{0,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}//Simplify;
AllComponentValues[Perturbation[EMT[\[Mu],-\[Nu]],0]//ToBasis[chart],(\[Rho][t[]]+p[t[]])u[\[Mu]]u[-\[Nu]]+p[t[]]metric[\[Mu],-\[Nu]]//SeparateMetric[metric]//ToBasis[chart]//ComponentArray//TraceBasisDummy//ToValues];
AllComponentValues[Perturbation[stress[\[Mu],-\[Nu]],0]//ToBasis[chart],{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}];
AllComponentValues[Perturbation[stress[\[Mu],-\[Nu]],1]//ToBasis[chart],stress1];
pertEMT1=(metric[\[Mu],-\[Nu]]\[Delta]p[t[],x[],y[],z[]]+Perturbation[metric[\[Mu],-\[Nu]],1]p[t[]]+Perturbation[u[\[Mu]],1]u[-\[Nu]](\[Rho][t[]]+p[t[]])+Perturbation[u[-\[Nu]],1]u[\[Mu]](\[Rho][t[]]+p[t[]])+u[\[Mu]]u[-\[Nu]](\[Delta]\[Rho][t[],x[],y[],z[]]+\[Delta]p[t[],x[],y[],z[]])+Perturbation[stress[\[Mu],-\[Nu]],1])//SeparateMetric[metric]//ExpandPerturbation//ToBasis[chart]//ComponentArray//TraceBasisDummy//ToValues//ToCanonical;
AllComponentValues[Perturbation[EMT[\[Mu],-\[Nu]],1]//ToBasis[chart],pertEMT1];


(* Define a scalar field phi and its potential V *)
DefTensor[phi[],M,PrintAs->"\[CurlyPhi]"];
DefTensorPerturbation[Pertphi[LI[order]],phi[],M,PrintAs->"\[Delta]\[CurlyPhi]"];
DefScalarFunction[V];

(* For gauge transformation of perturbations *)
DefScalarFunction[\[Xi]0];
DefScalarFunction[\[Xi]i];
DefTensor[gauge[\[Mu]],M];
DefTensorPerturbation[Pertgauge[LI[order],\[Mu]],gauge[\[Mu]],M];
PrintAs[Pertgauge]^="\[Xi]";
(*AllComponentValues[Perturbation[gauge[\[Mu]],0]//ToBasis[chart],{1,0,0,0}];*)
gaugeupP1=1/a[t[]] {a[t[]] \[Xi]0[t[],x[],y[],z[]],PDchart[{1,-chart}][\[Xi]i[t[],x[],y[],z[]]],PDchart[{2,-chart}][\[Xi]i[t[],x[],y[],z[]]],PDchart[{3,-chart}][\[Xi]i[t[],x[],y[],z[]]]};
AllComponentValues[Pertgauge[\[Xi]]//ToBasis[chart],gaugeupP1];


DefTensor[EinsteinEqleft[\[Mu],\[Nu]],M,Symmetric[{\[Mu],\[Nu]}]];
DefTensor[CartanEqleft[\[Sigma],-\[Mu],-\[Nu]],M,Antisymmetric[{-\[Mu],-\[Nu]}]];
DefTensor[EinsteinCompEqleftT[-\[Mu],-\[Nu]],M];
DefTensor[CartanCompEqleftT[\[Sigma],-\[Mu],-\[Nu]],M];


(* Pre-writing some time consuming terms *)
AllComponentValues[Perturbation[TorsionCD[#1,#2,#3],0]//ToBasis[chart],Perturbation[metric[#1,\[Tau]]metric[#2,\[Beta]]metric[#3,\[Gamma]]TorsionCD[-\[Tau],-\[Beta],-\[Gamma]],0]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues]&@@@Delete[Tuples[{{\[Sigma],-\[Sigma]},{\[Mu],-\[Mu]},{\[Nu],-\[Nu]}}],8];
Print["The components of TorsionCD have been written."];
Print[PGC121`bars];
AllComponentValues[Perturbation[Riccicd[#1,#2],0]//ToBasis[chart],Perturbation[metric[#1,\[Beta]]metric[#2,\[Gamma]]Riccicd[-\[Beta],-\[Gamma]],0]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues]&@@@Tuples[{{\[Mu],-\[Mu]},{\[Nu],-\[Nu]}}];
Print["The components of Riccicd have been written."];
Print[PGC121`bars];
AllComponentValues[ChristoffelcdPDchart[#1,#2,#3]//ToBasis[chart],(metric[#1,\[Tau]]metric[#2,\[Beta]]metric[#3,\[Gamma]]ChristoffelcdPDchart[-\[Tau],-\[Beta],-\[Gamma]])//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues]&@@@Tuples[{{\[Sigma],-\[Sigma]},{\[Mu],-\[Mu]},{\[Nu],-\[Nu]}}];
Print["The components of ChristoffelcdPDchart have been written."];
Print[PGC121`bars];
AllComponentValues[ChristoffelcdhPDchart[#1,#2,#3]//ToBasis[chart],(ProjectorToMetric[ChristoffelcdhPDchart[#1,#2,#3]//ChristoffelToMetric,metrich])//SeparateMetric[metric]//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues]&@@@Tuples[{{\[Sigma],-\[Sigma]},{\[Mu],-\[Mu]},{\[Nu],-\[Nu]}}];
Print["The components of ChristoffelcdhPDchart have been written."];
Print[PGC121`bars];
DefTensor[cdTorsionCD[\[Sigma],-\[Mu],-\[Nu],-\[Tau]],M,Antisymmetric[{-\[Mu],-\[Nu]}],PrintAs->"DT"];
DefTensor[cdcdTorsionCD[\[Sigma],-\[Mu],-\[Nu],-\[Tau],-\[Beta]],M,Antisymmetric[{-\[Mu],-\[Nu]}],PrintAs->"\!\(\*SuperscriptBox[\(D\), \(2\)]\)T"];
DefTensor[cdcdcdTorsionCD[\[Sigma],-\[Mu],-\[Nu],-\[Tau],-\[Beta],-\[Gamma]],M,Antisymmetric[{-\[Mu],-\[Nu]}],PrintAs->"\!\(\*SuperscriptBox[\(D\), \(3\)]\)T"];
DefTensor[cdcdRiemanncd[-\[Mu],-\[Nu],-\[Sigma],-\[Tau],-\[Beta],-\[Gamma]],M,RiemannSymmetric[{-\[Mu],-\[Nu],-\[Sigma],-\[Tau]}],PrintAs->"DDRie"];
DefTensor[cdRiccicd[-\[Sigma],-\[Mu],-\[Nu]],M,Symmetric[{-\[Sigma],-\[Mu]}],PrintAs->"DR"];
DefTensor[cdcdRiccicd[-\[Sigma],-\[Mu],-\[Nu],-\[Tau]],M,Symmetric[{-\[Sigma],-\[Mu]}],PrintAs->"\!\(\*SuperscriptBox[\(D\), \(2\)]\)R"];
DefTensor[cdRicciScalarcd[-\[Tau]],M,PrintAs->"DRs"];
DefTensor[cdcdRicciScalarcd[-\[Tau],-\[Sigma]],M,PrintAs->"DDRs"];
ruleDiffReplace={cd[\[Gamma]_]@cd[\[Beta]_]@cd[\[Tau]_]@TorsionCD[\[Sigma]_,\[Mu]_,\[Nu]_]->cdcdcdTorsionCD[\[Sigma],\[Mu],\[Nu],\[Tau],\[Beta],\[Gamma]],cd[\[Beta]_]@cd[\[Tau]_]@TorsionCD[\[Sigma]_,\[Mu]_,\[Nu]_]->cdcdTorsionCD[\[Sigma],\[Mu],\[Nu],\[Tau],\[Beta]],cd[\[Tau]_]@TorsionCD[\[Sigma]_,\[Mu]_,\[Nu]_]->cdTorsionCD[\[Sigma],\[Mu],\[Nu],\[Tau]],cd[\[Gamma]_]@cd[\[Beta]_]@Riemanncd[\[Mu]_,\[Nu]_,\[Sigma]_,\[Tau]_]->cdcdRiemanncd[\[Mu],\[Nu],\[Sigma],\[Tau],\[Beta],\[Gamma]],cd[\[Beta]_]@Riemanncd[\[Mu]_,\[Nu]_,\[Sigma]_,\[Tau]_]->cdRiemanncd[\[Mu],\[Nu],\[Sigma],\[Tau],\[Beta]],cd[\[Tau]_]@cd[\[Nu]_]@Riccicd[\[Sigma]_,\[Mu]_]->cdcdRiccicd[\[Sigma],\[Mu],\[Nu],\[Tau]],cd[\[Nu]_]@Riccicd[\[Sigma]_,\[Mu]_]->cdRiccicd[\[Sigma],\[Mu],\[Nu]],cd[\[Tau]_]@cd[\[Sigma]_]@RicciScalarcd[]->cdcdRicciScalarcd[\[Sigma],\[Tau]],cd[\[Sigma]_]@RicciScalarcd[]->cdRicciScalarcd[\[Sigma]]};
Module[{iter=0,length},
PrintTemporary["The components of cdTorsionCD are enumerating."];length=Length[Tuples[{{\[Sigma],-\[Sigma]},{\[Mu],-\[Mu]},{\[Nu],-\[Nu]},{\[Tau],-\[Tau]}}]];Monitor[Table[iter+=1;AllComponentValues[Perturbation[cdTorsionCD[#1,#2,#3,#4],0]//ToBasis[chart],Perturbation[cd[#4]@TorsionCD[#1,#2,#3],0]//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToValues]&@@Tuples[{{\[Sigma],-\[Sigma]},{\[Mu],-\[Mu]},{\[Nu],-\[Nu]},{\[Tau],-\[Tau]}}][[iter]],{iter,0,length-1}],ProgressIndicator[iter,{1,length}]];Print["The components of cdTorsionCD have been written."];];//AbsoluteTiming
Print[PGC121`bars];
Module[{},
PrintTemporary["The components of cdcdTorsionCD are enumerating."];AllComponentValues[Perturbation[cdcdTorsionCD[\[Sigma],-\[Mu],-\[Nu],-\[Tau],-\[Beta]],0]//ToBasis[chart],Perturbation[cd[-\[Beta]]@cdTorsionCD[\[Sigma],-\[Mu],-\[Nu],-\[Tau]],0]//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToValues];Print["The components of cdcdTorsionCD have been written."];];//AbsoluteTiming
Print[PGC121`bars];
Module[{},
PrintTemporary["The components of cdcdcdTorsionCD are enumerating."];AllComponentValues[Perturbation[cdcdcdTorsionCD[\[Sigma],-\[Mu],-\[Nu],-\[Tau],-\[Beta],-\[Gamma]],0]//ToBasis[chart],Perturbation[cd[-\[Gamma]]@cdcdTorsionCD[\[Sigma],-\[Mu],-\[Nu],-\[Tau],-\[Beta]],0]//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToValues];Print["The components of cdcdcdTorsionCD have been written."];];//AbsoluteTiming
Print[PGC121`bars];
Module[{iter=0,length},
PrintTemporary["The components of Riemanncd are enumerating."];length=Length[Tuples[{{\[Mu],-\[Mu]},{\[Nu],-\[Nu]},{\[Sigma],-\[Sigma]},{\[Tau],-\[Tau]}}]];Monitor[Table[iter+=1;AllComponentValues[Perturbation[Riemanncd[#1,#2,#3,#4],0]//ToBasis[chart],Perturbation[Riemanncd[#1,#2,#3,#4],0]//SeparateMetric[metric]//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToValues]&@@Tuples[{{\[Mu],-\[Mu]},{\[Nu],-\[Nu]},{\[Sigma],-\[Sigma]},{\[Tau],-\[Tau]}}][[iter]],{iter,0,length-1}],ProgressIndicator[iter,{1,length}]];Print["The components of Riemanncd have been written."];];//AbsoluteTiming
Print[PGC121`bars];
Module[{iter=0,length},
PrintTemporary["The components of cdcdRiemanncd are enumerating."];AllComponentValues[Perturbation[cdcdRiemanncd[-\[Mu],-\[Nu],-\[Sigma],-\[Tau],-\[Beta],-\[Gamma]],0]//ToBasis[chart],Perturbation[cd[-\[Gamma]]@cdRiemanncd[-\[Mu],-\[Nu],-\[Sigma],-\[Tau],-\[Beta]],0]//SeparateMetric[metric]//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToValues];Print["The components of cdcdRiemanncd have been written."];];//AbsoluteTiming
Print[PGC121`bars];
Module[{iter=0,length},
PrintTemporary["The components of cdRiccicd are enumerating."];length=Length[Tuples[{{\[Sigma],-\[Sigma]},{\[Mu],-\[Mu]},{\[Nu],-\[Nu]}}]];Monitor[Table[iter+=1;AllComponentValues[Perturbation[cdRiccicd[#1,#2,#3],0]//ToBasis[chart],Perturbation[cd[#3]@Riccicd[#1,#2],0]//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToValues]&@@Tuples[{{\[Sigma],-\[Sigma]},{\[Mu],-\[Mu]},{\[Nu],-\[Nu]}}][[iter]],{iter,0,length-1}],ProgressIndicator[iter,{1,length}]];Print["The components of cdRiccicd have been written."];];//AbsoluteTiming
Print[PGC121`bars];
Module[{iter=0,length},
PrintTemporary["The components of cdcdRiccicd are enumerating."];length=Length[Tuples[{{\[Sigma],-\[Sigma]},{\[Mu],-\[Mu]},{\[Nu],-\[Nu]},{\[Tau],-\[Tau]}}]];Monitor[Table[iter+=1;AllComponentValues[Perturbation[cdcdRiccicd[#1,#2,#3,#4],0]//ToBasis[chart],Perturbation[cd[#4]@cdRiccicd[#1,#2,#3],0]//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToValues]&@@Tuples[{{\[Sigma],-\[Sigma]},{\[Mu],-\[Mu]},{\[Nu],-\[Nu]},{\[Tau],-\[Tau]}}][[iter]],{iter,0,length-1}],ProgressIndicator[iter,{1,length}]];Print["The components of cdcdRiccicd have been written."];];//AbsoluteTiming
Print[PGC121`bars];
Module[{iter=0,length},
PrintTemporary["The components of cdRicciScalarcd are enumerating."];length=Length[Tuples[{{\[Sigma],-\[Sigma]}}]];Monitor[Table[iter+=1;AllComponentValues[Perturbation[cdRicciScalarcd[#1],0]//ToBasis[chart],Perturbation[cd[#1]@RicciScalarcd[],0]//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToValues]&@@Tuples[{{\[Sigma],-\[Sigma]}}][[iter]],{iter,0,length-1}],ProgressIndicator[iter,{1,length}]];Print["The components of cdRicciScalarcd have been written."];];//AbsoluteTiming
Print[PGC121`bars];
Module[{iter=0,length},
PrintTemporary["The components of cdcdRicciScalarcd are enumerating."];length=Length[Tuples[{{\[Sigma],-\[Sigma]},{\[Tau],-\[Tau]}}]];Monitor[Table[iter+=1;AllComponentValues[Perturbation[cdcdRicciScalarcd[#1,#2],0]//ToBasis[chart],Perturbation[cd[#2]@cdRicciScalarcd[#1],0]//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToValues]&@@Tuples[{{\[Sigma],-\[Sigma]},{\[Tau],-\[Tau]}}][[iter]],{iter,0,length-1}],ProgressIndicator[iter,{1,length}]];Print["The components of cdcdRicciScalarcd have been written."];];//AbsoluteTiming
Print[PGC121`bars];


(* the Lagrangian density and its perturbation *)
actionGeneralTo2th=\[Alpha] RicciScalarCD[]+a1 TorsionCD[\[Sigma],-\[Mu],-\[Nu]]TorsionCD[-\[Sigma],\[Mu],\[Nu]]+a2 TorsionCD[-\[Mu],-\[Nu],\[Sigma]]TorsionCD[\[Nu],\[Mu],-\[Sigma]]+a3 TorsionCD[\[Sigma],-\[Mu],-\[Sigma]]TorsionCD[\[Tau],\[Mu],-\[Tau]]+b1 RiemannCD[-\[Mu],-\[Nu],-\[Sigma],-\[Tau]]RiemannCD[\[Mu],\[Nu],\[Sigma],\[Tau]]+b2 RiemannCD[-\[Mu],-\[Nu],-\[Sigma],-\[Tau]]RiemannCD[\[Sigma],\[Tau],\[Mu],\[Nu]]+b3 RicciCD[-\[Mu],-\[Nu]]RicciCD[\[Mu],\[Nu]]+b4 RicciCD[-\[Mu],-\[Nu]]RicciCD[\[Nu],\[Mu]]+b5 RicciScalarCD[]^2+b6(epsilonmetric[-\[Mu],-\[Nu],-\[Sigma],-\[Tau]]RiemannCD[\[Mu],\[Nu],\[Sigma],\[Tau]])(epsilonmetric[\[Beta],\[Gamma],\[Xi],\[Zeta]]RiemannCD[-\[Beta],-\[Gamma],-\[Xi],-\[Zeta]])+c1 RiemannCD[-\[Mu],-\[Nu],-\[Sigma],-\[Tau]]RiemannCD[\[Mu],\[Sigma],\[Nu],\[Tau]]//ToCanonical;
coef={\[Alpha],a1,a2,a3,b1,b2,b3,b4,b5,b6,c1};
pert0Replace[Eqs_]:=Perturbation[Eqs,0]//.ruleDiffReplace;
pert1Replace[Eqs_]:=((Perturbation[Eqs,1]//ExpandPerturbation)//.ruleDeltaReplace//.ruleDiffReplace)//Expand//ToCanonical//ContractMetric;
parallel[pert_]:=If[MatchQ[pert,Y_ cdcdRiemanncd[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_,\[Gamma]_]]||MatchQ[pert,cdcdRiemanncd[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_,\[Gamma]_]]||MatchQ[pert,Y_ cdRiemanncd[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_]]||MatchQ[pert,cdRiemanncd[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_]]||MatchQ[pert,Y_ cdcdcdTorsionCD[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_,\[Gamma]_]]||MatchQ[pert,cdcdcdTorsionCD[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_,\[Gamma]_]]||MatchQ[pert,Y_ cdcdTorsionCD[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_]]||MatchQ[pert,cdcdTorsionCD[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_]]||MatchQ[pert,Y_ cdcdcdcddeltametric[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_,\[Gamma]_]]||MatchQ[pert,cdcdcdcddeltametric[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_,\[Gamma]_]]||MatchQ[pert,Y_ cdcdcddeltametric[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_]]||MatchQ[pert,cdcdcddeltametric[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_]]||MatchQ[pert,Y_ cdcdcddeltaTorsion[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_,\[Gamma]_]]||MatchQ[pert,cdcdcddeltaTorsion[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_,\[Gamma]_]]||MatchQ[pert,Y_ cdcddeltaTorsion[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_]]||MatchQ[pert,cdcddeltaTorsion[\[Sigma]_,\[Mu]_,\[Nu]_,\[Tau]_,\[Beta]_]],pert//SeparateMetric[metric]//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToValues//ToCanonical,pert//ToBasis[chart]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//ToValues//ToCanonical];
PrintTemporary["The field equations of action by term are reading."];
length=Length[coef];
time0=AbsoluteTime[];
Monitor[Do[iter+=1;
LagrangianDensity[]:=Sqrt[-Detmetric[]]Coefficient[actionGeneralTo2th,coef][[iter]];
PertLagrangian[]:=ToCanonical@ContractMetric@ExpandPerturbation@Perturbation@ToCanonical@(BreakChristoffel@NoScalar@ChangeCurvature[LagrangianDensity[],CD,cd]/.Which[$WhichVar==1,rule\[CapitalGamma]/.ruleK,$WhichVar==2,rule\[CapitalGamma]]);
(* Variations *)
If[!(FileExistsQ[FileNameJoin[{$PackageDirectory,StringTemplate["EinsteinEqleftAuxT`1`"][coef[[iter]]]}]]&&FileExistsQ[FileNameJoin[{$PackageDirectory,StringTemplate["CartanEqleftAuxT`1`"][coef[[iter]]]}]]),
EinsteinEqleft/:EinsteinEqleft[\[Mu]_,\[Nu]_]:=ToCanonical@ContractMetric@Symmetrize[metric[\[Mu],-\[Sigma]]metric[\[Nu],-\[Tau]]VarD[Perturbation[metric[-\[Sigma],-\[Tau]],1],cd][PertLagrangian[]]/Sqrt[-Detmetric[]]/.delta[-LI[1],LI[1]]->1//ToCanonical,{\[Mu],\[Nu]}];
CartanEqleft/:CartanEqleft[\[Sigma]_,\[Mu]_,\[Nu]_]:=ToCanonical@ContractMetric@Antisymmetrize[metric[\[Sigma],\[Tau]]metric[\[Mu],-\[Beta]]metric[\[Nu],-\[Gamma]]VarD[Perturbation[Which[$WhichVar==1,TorsionCD[\[Tau],-\[Beta],-\[Gamma]],$WhichVar==2,Contorsion[\[Tau],-\[Beta],-\[Gamma]]],1],cd][PertLagrangian[]]/Sqrt[-Detmetric[]]/.delta[-LI[1],LI[1]]->1//ToCanonical,Which[$WhichVar==1,{\[Mu],\[Nu]},$WhichVar==2,{\[Sigma],\[Nu]}]];
(* Symmetrization *)
EinsteinEqleftAuxT=EinsteinEqleft[-\[Mu],-\[Nu]]//ToCanonical;
CartanEqleftAuxT=CartanEqleft[\[Sigma],-\[Mu],-\[Nu]]//ToCanonical;

Put[EinsteinEqleftAuxT,FileNameJoin[{$PackageDirectory,StringTemplate["EinsteinEqleftAuxT`1`"][coef[[iter]]]}]];
Put[CartanEqleftAuxT,FileNameJoin[{$PackageDirectory,StringTemplate["CartanEqleftAuxT`1`"][coef[[iter]]]}]];,

EinsteinEqleftAuxT=Get[FileNameJoin[{$PackageDirectory,StringTemplate["EinsteinEqleftAuxT`1`"][coef[[iter]]]}]];
CartanEqleftAuxT=Get[FileNameJoin[{$PackageDirectory,StringTemplate["CartanEqleftAuxT`1`"][coef[[iter]]]}]];];

(* the components of Einstein and Cartan left parts *)
Off[General::argx];
(* the components of Equations of fields on the background*)
If[!FileExistsQ[FileNameJoin[{$PackageDirectory,StringTemplate["EinsteinCompEqleftT`1`"][coef[[iter]]]}]],
EinsteinCompEqleftT=Plus@@ParallelTable[parallel[#]&@pert0Replace[EinsteinEqleftAuxT][[ii]],{ii,1,Length@pert0Replace[EinsteinEqleftAuxT]},DistributedContexts->Automatic,Method->"FinestGrained"]//Simplification;
Put[EinsteinCompEqleftT,FileNameJoin[{$PackageDirectory,StringTemplate["EinsteinCompEqleftT`1`"][coef[[iter]]]}]];];

If[!FileExistsQ[FileNameJoin[{$PackageDirectory,StringTemplate["CartanCompEqleftT`1`"][coef[[iter]]]}]],
CartanCompEqleftT=If[iter==2||iter==3||iter==4,parallel[#]&@pert0Replace[CartanEqleftAuxT]//Simplification,Plus@@ParallelTable[parallel[#]&@pert0Replace[CartanEqleftAuxT][[ii]],{ii,1,Length@pert0Replace[CartanEqleftAuxT]},DistributedContexts->Automatic,Method->"FinestGrained"]//Simplification];
Put[CartanCompEqleftT,FileNameJoin[{$PackageDirectory,StringTemplate["CartanCompEqleftT`1`"][coef[[iter]]]}]];];

On[General::argx];,{iter,0,length-1}],ProgressIndicator[iter,{1,length}]];
time1=AbsoluteTime[];
Print[StringTemplate["The field equations of action by term have been written. Absolute timing: `1`seconds, `2`minutes, `3`hours."][time1-time0,(time1-time0)/60,(time1-time0)/3600]];


PrintTemporary["The Components of Lagrangian by term are reading."];
length=Length[coef];
time0=AbsoluteTime[];
Monitor[Do[iter+=1;
Off[General::argx];
If[!(FileExistsQ[FileNameJoin[{$PackageDirectory,StringTemplate["Lagrangian/LagrangianCompT`1`"][coef[[iter]]]}]]),
(* Symmetrization *)
LagrangianTerm=Coefficient[actionGeneralTo2th,coef][[iter]];
PreLagrangian=ToCanonical@ContractMetric@ToCanonical@(BreakChristoffel@NoScalar@ChangeCurvature[LagrangianTerm,CD,cd]/.rule\[CapitalGamma]/.ruleK);
LagrangianCompT=If[iter==2||iter==3||iter==4,parallel[#]&@pert0Replace[PreLagrangian]//Simplification,
Plus@@ParallelTable[parallel[#]&@pert0Replace[PreLagrangian][[ii]],{ii,1,Length@pert0Replace[PreLagrangian]},DistributedContexts->Automatic,Method->"FinestGrained"]//Simplification];

Put[LagrangianCompT,FileNameJoin[{$PackageDirectory,StringTemplate["Lagrangian/LagrangianCompT`1`"][coef[[iter]]]}]];];

On[General::argx];,{iter,0,length-1}],ProgressIndicator[iter,{1,length}]];
time1=AbsoluteTime[];
Print[StringTemplate["The Components of Lagrangian by term have been written. Absolute timing: `1`seconds, `2`minutes, `3`hours."][time1-time0,(time1-time0)/60,(time1-time0)/3600]];


EinsteinFieldEq[var_List]:=Plus@@Table[var[[iter]]Get[FileNameJoin[{$PackageDirectory,StringTemplate["EinsteinEqleftAuxT`1`"][var[[iter]]]}]],{iter,1,Length@var}]//ToCanonical;
CartanFieldEq[var_List]:=Plus@@Table[var[[iter]]Get[FileNameJoin[{$PackageDirectory,StringTemplate["CartanEqleftAuxT`1`"][var[[iter]]]}]],{iter,1,Length@var}]//ToCanonical;
EinsteinCompEqleftTvar[var_List]:=Plus@@Table[var[[iter]]Get[FileNameJoin[{$PackageDirectory,StringTemplate["EinsteinCompEqleftT`1`"][var[[iter]]]}]],{iter,1,Length@var}]//ToCanonical;
CartanCompEqleftTvar[var_List]:=Plus@@Table[var[[iter]]Get[FileNameJoin[{$PackageDirectory,StringTemplate["CartanCompEqleftT`1`"][var[[iter]]]}]],{iter,1,Length@var}]//ToCanonical;


(* the components equations *)
EinsteinCompEqDim[\[Mu]_Integer,\[Nu]_Integer,var_List]:=EinsteinCompEqleftTvar[var][[\[Nu]+1]][[\[Mu]+1]]==-kappa*(Perturbation[EMT[-\[Beta],-\[Gamma]]//SeparateMetric[metric],0]//ToBasis[chart]//TraceBasisDummy//ComponentArray//ToValues//Simplify)[[\[Nu]+1]][[\[Mu]+1]];
CartanCompEqDim[\[Mu]_Integer,\[Nu]_Integer,\[Rho]_Integer,var_List]:=CartanCompEqleftTvar[var][[\[Rho]+1]][[\[Nu]+1]][[\[Mu]+1]]==0;
EMTConservedCompEqDim[\[Mu]_Integer]:=(((ChangeCovD[CD[\[Beta]]@EMT[-\[Beta],-\[Gamma]],CD,cd]//BreakChristoffel)+TorsionCD[\[Tau],-\[Beta],-\[Tau]]EMT[\[Beta],-\[Gamma]]+TorsionCD[\[Beta],-\[Tau],-\[Gamma]]EMT[-\[Beta],\[Tau]])/.rule\[CapitalGamma]/.ruleK//ToCanonical//SeparateMetric[metric]//ToBasis[chart]//ToBasis[chart]//ComponentArray//TraceBasisDummy//ToValues//Simplify)[[\[Mu]+1]]==0;


LagrangianCompvar[var_List]:=Plus@@Table[var[[iter]]Get[FileNameJoin[{$PackageDirectory,StringTemplate["Lagrangian/LagrangianCompT`1`"][var[[iter]]]}]],{iter,1,Length@var}]//ToCanonical;


ruleH={PDchart[{0,-chart}][a[t[]]]->H[t[]] a[t[]],PDchart[{0,-chart}]@PDchart[{0,-chart}][a[t[]]]->(PDchart[{0,-chart}][H[t[]]]+H[t[]]^2) a[t[]],PDchart[{0,-chart}]@PDchart[{0,-chart}]@PDchart[{0,-chart}][a[t[]]]->(PDchart[{0,-chart}]@PDchart[{0,-chart}][H[t[]]]+3H[t[]] PDchart[{0,-chart}][H[t[]]]+H[t[]]^3) a[t[]],PDchart[{0,-chart}]@PDchart[{0,-chart}]@PDchart[{0,-chart}]@PDchart[{0,-chart}][a[t[]]]->(PDchart[{0,-chart}]@PDchart[{0,-chart}]@PDchart[{0,-chart}][H[t[]]]+H[t[]]^4+6H[t[]]^2 PDchart[{0,-chart}][H[t[]]]+4H[t[]] PDchart[{0,-chart}]@PDchart[{0,-chart}][H[t[]]]+3(PDchart[{0,-chart}][H[t[]]])^2) a[t[]](*,p[t[]]->\[Rho][t[]]w*)};


(*End[]*)


EndPackage[]
