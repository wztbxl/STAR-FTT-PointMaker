#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "TMath.h"
#include "TCernLib.h"

#include "StvSeedFinder.h"
#include "StvKalmanTrackFinder.h"

#include "StvToolkit.h"
#include "StMultiKeyMap.h"
#include "StarVMC/GeoTestMaker/StTGeoHelper.h"
#include "Stv/StvHit.h"
#include "StvUtil/StvNodePars.h"
#include "StvUtil/StvDebug.h"
#include "Stv/StvConst.h"
#include "Stv/StvDiver.h"
#include "Stv/StvHitter.h"
#include "Stv/StvFitter.h"
#include "Stv/StvTrackFitter.h"
#include "THelixTrack.h"
#include "Stv/StvDraw.h"
#include "Stv/StvStl.h"
#include "Stv/StvNode.h"
#include "Stv/StvTrack.h"
ClassImp(StvKalmanTrackFinder)
//_____________________________________________________________________________
StvKalmanTrackFinder::StvKalmanTrackFinder(const char *name):StvTrackFinder(name)
{
  memset(mBeg,0,mEnd-mBeg+1);
  mDive = new StvDiver("KalmanTrackFinder");
  mDive->Init();
  double rMax,zMin,zMax;
  StTGeoHelper::Inst()->GetHitShape()->Get(zMin,zMax,rMax);
  if (zMax < -zMin) zMax = -zMin;
  mDive->SetRZmax(rMax,zMax);
  mHitter = new StvHitter();
}  
//_____________________________________________________________________________
void StvKalmanTrackFinder::Clear(const char*)
{
 StvTrackFinder::Clear("");
}
//_____________________________________________________________________________
void StvKalmanTrackFinder::Reset()
{
}    

//_____________________________________________________________________________
int StvKalmanTrackFinder::FindTracks()
{
static StvToolkit *kit = StvToolkit::Inst();
static const StvKonst_st  *kons = StvConst::Inst();
  int nTrk = 0,nTrkTot=0,nAdded=0,nHits=0,nSeed=0,nSeedTot=0;
  StvSeedFinders *sfs = kit->SeedFinders();
  double aveRes=0,aveXi2=0;
  mCurrTrak = 0;

  for (int isf=0;isf<(int)sfs->size();isf++) { //Loop over seed finders
    StvSeedFinder* sf = (*sfs)[isf];
    int myMinHits = kons->mMinHits;
    if(sf->Again()) myMinHits = kons->mGoodHits;
    for (int repeat =0;repeat<5;repeat++) {//Repeat search the same seed finder 
      nTrk = 0;nSeed=0; sf->Again();
      while ((mSeedHelx = sf->NextSeed())) 
      {
	nSeed++;
						  if (sf->DoShow())  sf->Show();
	if (!mCurrTrak) mCurrTrak = kit->GetTrack();
	mCurrTrak->CutTail();	//Clean track from previous failure
	nAdded = FindTrack(0);
	if (!nAdded) 				continue;
//        StvNode *myNode = mCurrTrak->GetNode(StvTrack::kFirstPoint);
//        if (fabs(myNode->GetFP()._curv) <1./300) {
//  	  StvTrackFitter::Inst()->Helix(mCurrTrak,16|1);
//         }
	int ans = 0,fail=13;
    //		Refit track   
	do {
	  ans = Refit(1);
	  if (ans) 				break;
	  nAdded = FindTrack(1);
          if (nAdded<=0)			continue;;
    // few hits added. Refit track to beam again 
	  ans = Refit(0);
	  if (ans) 				break;
	} while((fail=0));		
	nHits = mCurrTrak->GetNHits();
	if (nHits < myMinHits)	fail+=100;		;
	if (fail) 	{//Track is failed, release hits & continue
	  mCurrTrak->CutTail();			continue;
        }
	StvNode *node = mCurrTrak->GetNode(StvTrack::kDcaPoint);
	if (node) node->UpdateDca();
        if (node && fabs(node->GetFP()._curv) <1./300) {
	  StvTrackFitter::Inst()->Helix(mCurrTrak,16|1);
        }

	kit->GetTracks().push_back(mCurrTrak);
	nTrk++;nTrkTot++;
	aveRes += mCurrTrak->GetRes();
	aveXi2 += mCurrTrak->GetXi2();
	mCurrTrak=0;
      }
      nSeedTot+=nSeed;
      Info("FindTracks:","SeedFinder(%s) Seeds=%d Tracks=%d ratio=%d\%\n"
          ,sf->GetName(),nSeed,nTrk,(100*nTrk)/(nSeed+1));
      
      if (!nTrk && myMinHits == kons->mMinHits) break;
      myMinHits = kons->mMinHits;
    }//End of repeat
  }//End of seed finders

  if (nTrkTot) {aveRes/=nTrkTot; aveXi2/=nTrkTot;}
  Info("FindTracks","tracks=%d aveRes = %g aveXi2=%g",nTrkTot,aveRes,aveXi2);
  return nTrkTot;
}
//_____________________________________________________________________________
int StvKalmanTrackFinder::FindTrack(int idir)
{

static int nCall=0; nCall++;
static const StvKonst_st *myConst = StvConst::Inst();
static       StvToolkit *kit      = StvToolkit::Inst();
static       StvFitter  *fitt     = StvFitter::Inst();

StvDebug::Break(nCall);
StvNodePars par[2];
StvFitErrs  err[2];
int mySkip=0,idive = 0,nNode=0,nHits=0;
double totLen=0;
StvNode *curNode=0,*preNode=0,*innNode=0,*outNode=0;
StvHitCount hitCount;
   
  if (mCurrTrak->empty()) {//Track empty, Backward tracking, to beam
    assert(!idir);
    double Hz = kit->GetHz(mSeedHelx->Pos());
    par[0].set(mSeedHelx,Hz); par[0].reverse();			//Set seed pars into par[0]
    err[0].Reset(Hz); 
  } else 	{//Forward or backward tracking
 
    curNode =(idir)? mCurrTrak->back(): mCurrTrak->front();
    par[0] = curNode->GetFP(); err[0] = curNode->GetFE(); 	//Set outer node pars into par[0]
  }

//  	Skip too big curvature
  if (fabs(par[0]._curv)>myConst->mMaxCurv)	return nHits;	

//  	skip P too small
  { double t = par[0]._tanl, pti = par[0]._ptin;	
    if ((t*t+1.)< myConst->mMinP2*pti*pti)	return nHits;
  }
  fitt->Set(par, err, par+1,err+1);
  fShowTrak.clear();
  mHitter->Reset();
StvFitDers derivFit;
  mDive->Reset();
  mDive->Set(par+0,err+0,idir);
  mDive->Set(par+1,err+1,&derivFit);	//Output of diving in par[1]


  while(idive==kDiveHits || idive==0) {

    do {//Stop tracking?
      idive = 99;
      if (!nNode)		continue;
      mySkip = hitCount.Skip();
      if (!mySkip) 		continue;
      if (idir)   		break;
      if (hitCount.Reject()) 	break;
      mDive->SetSkip();
    } while ((idive=0));
    if (idive) 			break;

//+++++++++++++++++++++++++++++++++++++

    idive = mDive->Dive();

//+++++++++++++++++++++++++++++++++++++
    if (mySkip && !idive) {
      Warning("FindTrack","Strange case mySkip!=0 and iDive==0");
      break;
    }
    totLen+=mDive->GetLength();
    nNode++;		// assert(nNode<200);
    if (nNode>200) { //Something very wrong
      Error("FindTrack","Too many nodes =200 Skip track");
      return 0;
    }
    if (idive >= kDiveBreak) 			break;
    par[0]=par[1]; err[0]=err[1];			//pars again in par[0]
		// Stop tracking when too big Z or Rxy
    if (fabs(par[0]._z)  > myConst->mZMax  ) 	break;
    if (par[0].getRxy()  > myConst->mRxyMax) 	break;

    		
    if (err[0].Check("AfterDive"))		break;
//    assert(idive || !err[0].Check("AfterDive"));
    const StvHits *localHits = 0; 
    if (idive== kDiveHits) {
static float gate[2]={myConst->mMaxWindow,myConst->mMaxWindow};   
      localHits = mHitter->GetHits(par,gate); 
    }

if (DoShow()) {
//  double r = par[0].getRxy();
    printf("%3d Len=%g XY=%g %g Path=%s\n",nNode,totLen,par[0]._x,par[0]._y
          ,StTGeoHelper::Inst()->GetPath());
    fShowTrak+=&par[0]._x;
}//EndDoShow

//	Create and add nodemyTrak
    preNode = curNode;
    curNode = kit->GetNode();      
    if (!idir)  {mCurrTrak->push_front(curNode);innNode=curNode;outNode=preNode;}
    else        {mCurrTrak->push_back (curNode);innNode=preNode;outNode=curNode;}

//    assert(!idive || !par[0].check("FindTrack.1"));
    curNode->mLen = (!idir)? totLen:-totLen;
		// Set prediction
    curNode->SetPre(par[0],err[0],0);
    assert(!preNode || innNode->GetFP().getRxy()<outNode->GetFP().getRxy());
    innNode->SetDer(derivFit,idir);
    innNode->SetELoss(mDive->GetELossData(),idir);

    if (idive==kDiveDca) {
      curNode->SetType(StvNode::kDcaNode);
      double testDca = TCL::vdot(&par[0]._cosCA,par[0].P,2);
      assert(fabs(testDca)<1e-4);
      continue;
    }

    if (!localHits)	 continue;	//Never hits in curNode 
    curNode->SetHitPlane(mHitter->GetHitPlane());
    if (!localHits->size()) {//No hits in curNode
      hitCount.AddNit(); continue;
    } 
    fitt->Prep();
    double minXi2 = myConst->mXi2Hit; 
    StvHit *minHit=0; int minIdx = -1;
    for (int ihit=0;ihit<(int)localHits->size();ihit++) {
      StvHit *hit = (*localHits)[ihit];
      double myXi2 = fitt->Xi2(hit);
      if (myXi2 > minXi2) continue;
      minXi2 = myXi2; minHit = hit; minIdx = ihit;
    }
    curNode->SetHit(minHit); 
    if (! minHit) {
      hitCount.AddNit(); 
    } else {
      hitCount.AddHit();

      (*((StvHits*)localHits))[minIdx]=0;
      double myXi2 = fitt->Xi2(minHit);
      assert(fabs(minXi2-myXi2)<1e-5);
		//fitted pars again in par[1]
      int iuerr = fitt->Update();if (iuerr){}; 
      nHits++;
      if (nHits<=3) par[1]=par[0]; //Fit is not reliable yet
      if (par[1].check()) {	//Ugly parameters
        hitCount.AddNit(); nHits--; curNode->SetHit(0);}
      else {
        hitCount.AddHit();
        curNode->SetXi2(myXi2,0);
        curNode->SetHE(fitt->GetHitErrs());
        curNode->SetFit(par[1],err[1],0);
      }

		// par[0] again keeps the latest version og pars
      par[0]=par[1]; err[0]=err[1];
    }



  } // End Dive&Fitter loop 

  mCurrTrak->SetTypeEnd(mySkip);
  if (!idir) {
    double eff = hitCount.Eff();
    if (hitCount.Reject()) {
      StvDebug::Count("BadEff",eff);
      if (hitCount.nContNits)StvDebug::Count("BadCNits",hitCount.nContNits);
      if (hitCount.nTotNits )StvDebug::Count("BadTNits",hitCount.nTotNits);
      mCurrTrak->ReleaseHits(); mCurrTrak->unset();
      kit->FreeTrack(mCurrTrak);mCurrTrak=0; return 0; }
    if (hitCount.nContNits)StvDebug::Count("GooCNits",hitCount.nContNits);
    if (hitCount.nTotNits )StvDebug::Count("GooTNits",hitCount.nTotNits);
    StvDebug::Count("GoodEff",eff);
  }



  return nHits;

}
//_____________________________________________________________________________
int StvKalmanTrackFinder::FindPrimaries(const StvHits &vtxs)	
{
static const StvKonst_st *myConst =   StvConst::Inst();
static     StvToolkit *kit     = StvToolkit::Inst();
static     StvTrackFitter *tkf = StvTrackFitter::Inst();

  StvTracks &traks = kit->GetTracks();
  int goodCount= 0, plus=0, minus=0;
  int nVertex = vtxs.size();  
  if (!nVertex) return 0;
  int nTracks = 0;
  for (StvTrackIter it=traks.begin(); it!=traks.end() ;++it) {
    StvTrack *track = *it;  nTracks++;
    double dca00 = track->ToBeam();
    if (dca00 > myConst->mDca2dZeroXY) {
      if (dca00 >1e11) StvDebug::Count("PrimNoDcaRej",    0);
      else             StvDebug::Count("PrimDca00Rej",dca00);
      continue;
    }
    int bestVertex=-1; double bestXi2 = myConst->mXi2Vtx;
    for (int iVertex=0;iVertex<nVertex;iVertex++) {
      StvHit *vertex = vtxs[iVertex];
      if (tkf->Fit(track,vertex,0)) 		continue;
      double Xi2 = tkf->GetXi2();
      if (Xi2>=bestXi2) 			continue;
// 		Found better Xi2
      bestXi2 = Xi2; bestVertex=iVertex;
    }//End vertex loop
    
    if(bestVertex<0) 				continue;
    StvDebug::Count("PrimXi2Acc",bestXi2);
    StvNode *node = kit->GetNode();
    StvHit *hit = vtxs[bestVertex];
    hit->addCount();
    tkf->Fit(track,hit,node);
    track->push_front(node);
    track->SetPrimary(bestVertex+1);
    node->SetType(StvNode::kPrimNode);    
    node->SetHit(hit);    
    node->SetXi2(bestXi2,0);
    goodCount++;
    if (track->GetCharge()>0) { plus++; } else { minus++; }

  }//End track loop 
  return goodCount;
}
//_____________________________________________________________________________
int StvKalmanTrackFinder::Refit(int idir)
{
static int nCall=0;nCall++;
static StvTrackFitter *tkf = StvTrackFitter::Inst();
  StvNode *node = 0;
  int ans=0;
  int lane = 1;
  for (int refIt=0; refIt<5; refIt++)  {

    ans = tkf->Refit(mCurrTrak,idir,lane,1);
    if (!ans) 		break;	//SUCCESS
    int nHits = mCurrTrak->GetNHits();
    if (nHits<3) 	return 10;

    printf("AfterRefit: worst Xi2 %g(%p)\n"
          ,tkf->GetWorstXi2(),tkf->GetWorstNode());

    tkf->Helix(mCurrTrak,1|2);
    node = tkf->GetWorstNode();

    printf("AfterHelix: worst Xi2 %g(%p)\n"
          ,tkf->GetWorstXi2(),tkf->GetWorstNode());

    assert(node);
    node->SetHit(0);
    ans = tkf->Refit(mCurrTrak,1-idir,1-lane,1);
  }
  node = mCurrTrak->GetNode(StvTrack::kDcaPoint);
  if (!node) return ans;
  node->UpdateDca();
  return ans;
}

