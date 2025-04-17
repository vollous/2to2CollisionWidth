(* ::Package:: *)

(* ::Input:: *)
(*$LoadAddOns={"FeynArts"};*)
(*<<FeynCalc`*)
(*$FAVerbose=0;*)


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
(*ampSquared[0]= (amp[0] (ComplexConjugate[amp[0]]))//FeynAmpDenominatorExplicit//SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#]&//DoPolarizationSums[#,p2,0]&//DoPolarizationSums[#,k1,0]&//DiracSimplify//TrickMandelstam[#,{s,t,u,2 SMP["m_t"]^2+ms^2}]&//FullSimplify*)


(* ::Input:: *)
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +mg^2+ms^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->e}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


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
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +mg^2+ms^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->e}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


(* ::Subsection:: *)
(*Subscript[t, L ] Z\[RightArrow] Subscript[t, R]Z*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[2]}->{V[2],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]];
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
(*ampSquared[1]=ampSquared[0]//.{u->2SMP["m_t"]^2 +mg^2+ms^2-s - t,SUNN->3,SMP["m_t"]->mt,SMP["m_W"]->mw,SMP["sin_W"]->sw,SMP["g_s"]->gs,SMP["e"]->e}*)


ampSquared[2]=(Asymptotic[Asymptotic[((ampSquared[1]//Numerator)/.{ms->0}),mg->0],mt->0]/.{mt->mtpole}) / (ampSquared[1]//Denominator)//Simplify


 


StringReplace[ToString[ampSquared[2]//CForm],{"Power"->"pow","mtpole"->"mt_pole","mt,2"->"mtinf,2","sw"->"sW","e"->"el","mw"->"mW"}]


(* ::Subsection:: *)
(*Subscript[t, L ] W\[RightArrow] Subscript[t, R]W*)


diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],V[3]}->{V[3],F[3,{3}]},InsertionLevel->{Classes},Model->"SMQCD"][[{2}]]
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
