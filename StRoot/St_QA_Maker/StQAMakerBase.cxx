// $Id: StQAMakerBase.cxx,v 2.7 2001/05/23 00:14:53 lansdell Exp $ 
// $Log: StQAMakerBase.cxx,v $
// Revision 2.7  2001/05/23 00:14:53  lansdell
// more changes for qa_shift histograms
//
// Revision 2.6  2001/05/16 20:57:03  lansdell
// new histograms added for qa_shift printlist; some histogram ranges changed; StMcEvent now used in StEventQA
//
// Revision 2.5  2001/04/28 22:05:13  genevb
// Added EMC histograms
//
// Revision 2.4  2001/04/26 16:34:40  genevb
// Fixed some histogram ranges
//
// Revision 2.3  2000/08/25 22:06:15  genevb
// Raised high mult bin to 2500 global tracks
//
// Revision 2.2  2000/08/25 20:29:34  lansdell
// year1 print list changed slightly; cosmetic improvement of some event summary histograms
//
// Revision 2.1  2000/08/25 16:04:10  genevb
// Introduction of files
//
// 
///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  StQAMakerBase base class for QA Histogram Makers                     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "StMessMgr.h"
#include "StQAMakerBase.h"
#include "StQABookHist.h"
#include "QAH.h"
#include "TList.h"

ClassImp(StQAMakerBase)

//_____________________________________________________________________________
StQAMakerBase::StQAMakerBase(const char *name, const char *title, const char* type) :
 StMaker(name,title), QAMakerType(type) {

  hists = 0;
  histsList = new TList();
  histsSet = 0;
  multClass = 3;

//  - Set all the histogram booking constants

  ntrk   = 50;
  nmnpt  = 50;
  nmneta = 50;
  nxyz   = 50;

//  - Zero all histogram pointers
  mNullPrimVtx = 0; // histogram for number of events without primary vertex
  mMultClass = 0;   // histogram for number of events in mult classes

// for method MakeEvSum - from table software monitor
  m_glb_trk_chg=0;          //! all charge east/west, tpc
  m_glb_trk_chgF=0;         //! all charge east/west, ftpc

}
//_____________________________________________________________________________
StQAMakerBase::~StQAMakerBase() {
  delete histsList;
}
//_____________________________________________________________________________
Int_t StQAMakerBase::Init(){
// Histogram booking must wait until first event Make() to determine event type
  firstEvent = kTRUE;
  return StMaker::Init();
}
//_____________________________________________________________________________
Int_t StQAMakerBase::Make(){

  gMessMgr->Info(" In StQAMakerBase::Make()");

  if (firstEvent) BookHist();
  // See BookHist() for definitions of histsSet values,
  // which should be set during Make() of the derived QA Maker class
  switch (histsSet) {
    case (0) : {
      multClass = 0;
      hists = (StQABookHist*) histsList->At(0);
      break; }
    case (1) : {
      if (multiplicity < 50) multClass = 0;
      else if (multiplicity < 500) multClass = 1;
      else if (multiplicity < 2500) multClass = 2;
      else multClass = 3;

      mMultClass->Fill(multClass);
      if (!multClass) return kStOk;
      hists = (StQABookHist*) histsList->At((Int_t) (--multClass));
      break; }
    default  : {}
  }

  // Call methods to fill histograms

  // histograms from table event_summary
  MakeHistEvSum();
  // histograms from table globtrk
  MakeHistGlob();
  // histograms from table primtrk - must be done after global tracks
  MakeHistPrim();
  // histograms from table primtrk & dst_dedx
  MakeHistPID();
  // histograms from table dst_dedx
  MakeHistDE();
  // histograms from table point
  MakeHistPoint();
  // histograms from table g2t_rch_hit
  MakeHistRich();
  // histograms from table dst_vertex,dst_v0_vertex,dst_xi_vertex,dst_kinkVertex
  MakeHistVertex();
  // histograms from EMC in StEvent
  MakeHistEMC();
  // histograms from geant and reco tables 
  if (histsSet==0) MakeHistEval();

  return kStOk;
}
//_____________________________________________________________________________
void StQAMakerBase::NewQABookHist(const char* prefix) {  

  gMessMgr->Info() <<
    "StQAMakerBase: booking histograms with prefix: " << prefix << endm;
  hists = new StQABookHist(prefix);
  histsList->Add(hists);
}
//_____________________________________________________________________________
TH2F* StQAMakerBase::MH1F(const Text_t* name, const Text_t* title,
                          Int_t nbinsx, Axis_t xlow, Axis_t xup) {
  TH2F* h = QAH::MH1F(name,title,nbinsx,xlow,xup,(Int_t) multClass);
  if (multClass>1) {
    h->Rebin(0,"low mult");
    h->Rebin(1,"mid mult");
    h->Rebin(2,"high mult");
  }
  return h;
}
//_____________________________________________________________________________
void StQAMakerBase::BookHist() {  
// book histograms

  firstEvent = kFALSE;
  TString temp;

  switch (histsSet) {

    // Generic data (e.g. Monte Carlo) with just one multiplicity class
    case (0) : {
      NewQABookHist(QAMakerType.Data());
      break; }

    // Real data with three multiplicity classes (low, medium, high)
    case (1) : {
      (temp = QAMakerType) += "LM";
      NewQABookHist(temp.Data());
      (temp = QAMakerType) += "MM";
      NewQABookHist(temp.Data());
      (temp = QAMakerType) += "HM";
      NewQABookHist(temp.Data());
      break; }

    default  : {}
  }
  
  multClass = histsList->GetSize();
  QAH::maker = (StMaker*) (this);
  QAH::preString = QAMakerType;

  BookHistGeneral();
  BookHistEvSum();
  for (Int_t i=0; i<multClass; i++)
    ((StQABookHist*) (histsList->At(i)))->BookHist(histsSet);
}
//_____________________________________________________________________________
void StQAMakerBase::BookHistGeneral(){  

  mNullPrimVtx = QAH::H1F("QaNullPrimVtx","event primary vertex check",40,-2,2);
  mNullPrimVtx->SetXTitle("has primary vertex? (yes = 1, no = -1)");
  mNullPrimVtx->SetYTitle("# of events");

  if (histsSet != 0) {
    mMultClass = QAH::H1F("QaMultClass","event multiplicity class",5,-0.5,4.5);
    mMultClass->SetXTitle("mult class (0=?/MC, 1=LM, 2=MM, 3=HM)");
    mMultClass->SetYTitle("# of events");
  }
}
//_____________________________________________________________________________
void StQAMakerBase::BookHistEvSum(){  

// for method MakeEvSum - from software monitor

  m_glb_trk_chg = MH1F("QaEvsumTotChg","softmon: all charge east/west,tpc",60,0,3);
  m_glb_trk_chgF = MH1F("QaEvsumTotChgF","softmon: all charge east/west,ftpc",60,0,3);
}
//_____________________________________________________________________________
