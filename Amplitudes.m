(* ::Package:: *)

(* ::Input:: *)
(*$LoadAddOns={"FeynArts"};*)
(*<<FeynCalc`*)
(*$FAVerbose=0;*)


(* ::Input:: *)
(*diags=InsertFields[CreateTopologies[0,2->2],{F[3,{3}],-F[3,{3}]}->{V[5],V[5]},InsertionLevel->{Classes},Model->"SMQCD"];Paint[diags,ColumnsXRows->{2,1},Numbering->Simple,SheetHeader->None,ImageSize->{512,256}];*)


(* ::Input:: *)
(*amp[0]=FCFAConvert[CreateFeynAmp[diags],IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,TransversePolarizationVectors->{k1,k2},List->False,SMP->True,Contract->True,DropSumOver->True]*)


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
