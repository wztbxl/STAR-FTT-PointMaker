// $Id: StKFVertexMaker.h,v 2.2 2013/04/08 19:21:41 fisyak Exp $

#ifndef STAR_StKFVertexMaker
#define STAR_StKFVertexMaker

/*!
 *                                                                     
 * \class  StKFVertexMaker
 * \author fisyak
 * \date   2012/04/18
 * \brief  virtual base class for Maker
 *
 */                                                                      

#ifndef StMaker_H
#include "StMaker.h"
#endif
#include "TObjArray.h"
#include "TSpectrum.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"
#include "StEnumerations.h"
#include "TCanvas.h"
class StEvent;
class StPrimaryVertex;
class StEvent;
class StGlobalTrack;
class StDcaGeometry;
class StTrack;
class KFParticle;
class StKFVerticesCollection;
class StiHit;
class StiKalmanTrackNode;
class StiKalmanTrack;
class StKFVertexMaker : public StMaker {
 public: 
  StKFVertexMaker(const char *name="KFVertex");
  virtual       ~StKFVertexMaker();
  virtual Int_t Init();
  virtual Int_t Make();
  // virtual Int_t FinishRun(Int_t runumber){return 0;}; // Overload empty StMaker::FinishRun 
  void Clear(Option_t *option="");
  void Fit();
  TH1F *VtxM() {return fVtxM;}
  void SetZwindow(Double_t z = 2) {fzWindow = z;}
  void SetDefaultTempLog(Double_t tLog = 2) {fTempLog = tLog;}
  static Double_t AnnelingFcn(Double_t TInv=1);
  TH1 *Vtx() {return fVtx;}
  StKFVerticesCollection* Vertices() {return fcVertices;}
  TObjArray &Particles() {return *fParticles;}
  KFParticle *AddTrackAt(const StGlobalTrack *gTrack);
  KFParticle *AddTrackAt(const StDcaGeometry *dca,Int_t kg);
  KFParticle *AddTrackAt(const StiKalmanTrackNode *tNode,Int_t kg);
  KFParticle *AddBeamTrack(); 
  void calculateRank(StPrimaryVertex *primV);
  void SetCanvas(TCanvas *c1) {fc1 = c1;}
  TCanvas *Canvas() {return fc1;}
  TH1F *GetVtxs(Int_t pass = 0) {return fVtxs[pass];}
  TH1K *GetVtxKs(Int_t pass = 0) {return fVtxKs[pass];}
  TH1F *GetVtxM() {return fVtxM;}
  Int_t MakeParticles();
  void  UpdateParticleAtVertex(StiKalmanTrack */* kTrack */, KFParticle */* particle */);
  Bool_t MakeV0(StKFVertex *V);
  StTrack  *FitTrack2Vertex(StKFVertex *V, const StKFTrack*   track);
 private:
  TObjArray *fParticles; // KF particles
  Int_t fNzBins;
  Int_t fNPasses;
  TSpectrum *fSpectrum;
  Double_t fzWindow;
  TH1F  *fVtxM;
  StKFVerticesCollection **fVerticesPass;
  static StKFVerticesCollection *fcVertices;  // current vertex collection
  Double_t fTempLog;
  ROOT::Math::GSLMinimizer1D *fminBrent;
  ROOT::Math::Functor1D      *func;
  TH1F *fVtxs[2];
  TH1  *fVtx;
  TH1K *fVtxKs[2];
  Bool_t mBeamLine;
  StPrimaryVertexOrder     mVertexOrderMethod; // will default to 0 i.e. orderByNumberOfDaughters
  TCanvas                 *fc1;
  StEvent *pEvent;
  /// Displayed on session exit, leave it as-is please ...
  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $ $Id: StKFVertexMaker.h,v 2.2 2013/04/08 19:21:41 fisyak Exp $ built "__DATE__" "__TIME__ ; 
    return cvs;
  }

  ClassDef(StKFVertexMaker,0)   //StAF chain virtual base class for Makers
};
#endif
// $Log: StKFVertexMaker.h,v $
// Revision 2.2  2013/04/08 19:21:41  fisyak
// Adjust for new KFParticle
//
// Revision 2.1  2012/05/07 14:56:14  fisyak
// Add StKFVertexMaker
//
