(* ::Package:: *)

(* ::Input:: *)
(*$LoadAddOns={"FeynArts"};*)
(*<<FeynCalc`*)
(*$FAVerbose=0;*)

commonpart=\!\(TraditionalForm\`
\*FractionBox[\(\(\ \)\(
\*SuperscriptBox[\(mtpole\), \(2\)]\ \((s + t)\)\)\), \(\(\ \)
\*SuperscriptBox[\((
\*SuperscriptBox[\(mt\), \(2\)] - t)\), \(2\)]\(\ \)\)]\);

commonpartMassive=\!\(TraditionalForm\`
\*FractionBox[\(\(\ \)\(
\*SuperscriptBox[\(mtpole\), \(2\)]\ \((s + t)\) t^2\)\), \(\(\ \)
\*SuperscriptBox[\((
\*SuperscriptBox[\(mt\), \(2\)] - t)\), \(2\)]\(\ \)\)]\);
expvalues={sw->Sqrt[0.223],SMP["cos_W"]->Sqrt[1-sw^2],EE->0.332,gs->Sqrt[4 \[Pi] alphaS],alphaS->0.12};


(* ::Subsection:: *)
(*Subscript[t, L ]g \[RightArrow] Subscript[t, R] H*)


(* ::Input:: *)
(*diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[5]}->{S[1],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];*)


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],0,ms,SMP["m_t"]];*)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2,0]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+ms^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +mg^2+ms^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->e}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)


 ampSquared[2]/commonpart
 %//.expvalues//FullSimplify


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


(* ::Subsection:: *)
(*Subscript[t, L ] g \[RightArrow] Subscript[t, R] g*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[5]}->{V[5],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];
Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2,k1},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],0,0,SMP["m_t"]];*)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2,0]&//DoPolarizationSums[#,k1,0]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+2mg^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +2 mg^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->EE}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


 ampSquared[2]/commonpart
 %//.expvalues


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


(* ::Subsection:: *)
(*Subscript[t, L ] g \[RightArrow] Subscript[t, R] \[Gamma]*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[5]}->{V[1],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];
Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2,k1},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],0,0,SMP["m_t"]];*)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2,0]&//DoPolarizationSums[#,k1,0]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+ms^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +mg^2+ms^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->EE}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


 ampSquared[2]/commonpart
 %//.expvalues


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


(* ::Subsection:: *)
(*Subscript[t, L ] \[Gamma] \[RightArrow] Subscript[t, R]\[Gamma]*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[1]}->{V[1],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];
Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2,k1},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],0,0,SMP["m_t"]];*)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2,0]&//DoPolarizationSums[#,k1,0]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+ms^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +mg^2+ms^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->EE}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


 ampSquared[2]/commonpart
 %//.{sw->Sqrt[0.223],SMP["cos_W"]->Sqrt[1-sw^2],EE->0.332}


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


(* ::Subsection:: *)
(*Subscript[t, L ] Z\[RightArrow] Subscript[t, R]Z*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[2]}->{V[2],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];
Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2,k1},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SP[k1,k1]=mz^2;*)
(*SP[p2,p2]=mz^2;*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],mz,mz,SMP["m_t"]];*)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2]&//DoPolarizationSums[#,k1]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+mz^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +mg^2+ms^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->EE}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mz->0],mt->0]/.{mt->mtpole,mg->0}) / (ampSquared[1]//Denominator)//Simplify


 ampSquared[2]/commonpart
 %//.{sw->Sqrt[0.223],SMP["cos_W"]->Sqrt[1-sw^2],EE->0.332}


(* ::Subsection:: *)
(*Subscript[t, L ] h\[RightArrow] Subscript[t, R]h*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],S[1]}->{S[1],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];
Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2,k1},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],0,0,SMP["m_t"]];*)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2,0]&//DoPolarizationSums[#,k1,0]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+ms^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +mg^2+ms^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->EE}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


 ampSquared[2]/commonpart
 %//.expvalues


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


(* ::Subsection:: *)
(*Subscript[t, L ] W\[RightArrow] Subscript[t, R]W*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[3]}->{F[3,{3}],V[3]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]]
Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2,k2},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],SMP["m_W"],SMP["m_t"],SMP["m_W"]];*)
(**)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2]&//DoPolarizationSums[#,k2]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+2mw^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +2mw^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_Z"]->mz,SMP["m_H"]->mh,SMP["m_W"]->mw,mw->Sqrt[1-sw^2] mz,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->e}//Expand//Simplify*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{mz->0}),mt->0],mt->2]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


(Asymptotic[Asymptotic[Asymptotic[((ampSquared[1]//Numerator)),mz->0],mt->0],mh->0]//Simplify)/ (ampSquared[1]//Denominator)//Simplify


(* ::Subsection:: *)
(*Subscript[t, L ] g \[RightArrow] Subscript[t, R] \[Gamma]*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[5]}->{V[1],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];
Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2,k1},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],0,0,SMP["m_t"]];*)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2,0]&//DoPolarizationSums[#,k1,0]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+ms^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +mg^2+ms^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->e}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


 


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


alphaS  = 0.12





64 Sqrt[4 * Pi * alphaS]^2/3


ee =0.332;


Sqrt[4 * Pi * alphaS]^2


64/3//N


64/27 * 0.332^4


((4 * ee^4 (3-4 sw^2)^2)/(27 Cos[theta]^4))//.{sw->Sqrt[0.223],theta->ArcSin[sw]}


3*ee^4/(4 sw^4)/.sw->Sqrt[0.223]


(64 e^2 gs^2)/9/. {e->ee,gs->Sqrt[4 \[Pi] alphaS]}


32/1.18


(* ::Title:: *)
(*Massive Cases*)


(* ::Subsection:: *)
(*Subscript[t, L ] g \[RightArrow] Subscript[t, R] g*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[5]}->{V[5],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];
Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2,k1},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],mg,mg,SMP["m_t"]];*)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2]&//DoPolarizationSums[#,k1]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+2mg^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +2 mg^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->EE}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


 ampSquared[2]/commonpartMassive
 %//.expvalues


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


(* ::Subsection:: *)
(*Subscript[t, L ] g \[RightArrow] Subscript[t, R] \[Gamma]*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[5]}->{V[1],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];
Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2,k1},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],mg,mphoton,SMP["m_t"]];*)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2]&//DoPolarizationSums[#,k1]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+mphoton^2+mg^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +mphoton^2+ mg^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->EE}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{mphoton->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


 ampSquared[2]/commonpartMassive
 %//.expvalues


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


(* ::Subsection:: *)
(*Subscript[t, L ] g \[RightArrow] Subscript[t, R] \[Gamma] (massless)*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[5]}->{V[1],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];
Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2,k1},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],mg,0,SMP["m_t"]];*)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2]&//DoPolarizationSums[#,k1,0]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+mg^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 + mg^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->EE}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{mphoton->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


 ampSquared[2]/commonpartMassive
 %//.expvalues


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


(* ::Subsection:: *)
(*Subscript[t, L ] \[Gamma] \[RightArrow] Subscript[t, R] \[Gamma]*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[1]}->{V[1],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];
Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{p2,k1},List->False,SMP->True,Contract->True,DropSumOver->True]/.{Spinor[Momentum[k2],SMP["m_t"],1]->Spinor[Momentum[k2],SMP["m_t"],1] . DiracGamma[7],Spinor[Momentum[p1],SMP["m_t"],1]->DiracGamma[7] . Spinor[Momentum[p1],SMP["m_t"],1]}*)


(* ::Input:: *)
(*FCClearScalarProducts[];*)
(*SetMandelstam[s,t,u,p1,p2,-k1,-k2,SMP["m_t"],mphoton,mphoton,SMP["m_t"]];*)


(* ::Input:: *)
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2]&//DoPolarizationSums[#,k1]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+2mphoton^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +2 mphoton^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->EE}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mphoton->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


 ampSquared[2]/commonpartMassive
 %//.expvalues


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]
