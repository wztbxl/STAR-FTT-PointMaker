/***************************************************************************
 *
 * $Id: StTofSimMaker.cxx,v 1.1 2001/09/28 19:11:11 llope Exp $
 *
 * Author: Frank Geurts
 ***************************************************************************
 *
 * Description: StTofSimMaker class for TOFp Simulations
 *
 ***************************************************************************
 *
 * $Log: StTofSimMaker.cxx,v $
 * Revision 1.1  2001/09/28 19:11:11  llope
 * first version
 *
 **************************************************************************/
#include <iostream.h>
#include <iomanip.h>    //fg only now needed for setw()
#include "StTofSimMaker.h"

// SCL
#include "StThreeVectorD.hh"
#include "Random.h"
#include "RanluxEngine.h"
#include "RandGauss.h"

#include "StTofUtil/StTofCalibration.h"
#include "StTofUtil/StTofSimParam.h"
#include "StTofUtil/StTofGeometry.h"
#include "StEventTypes.h"

// g2t tables
#include "tables/St_g2t_ctf_hit_Table.h"
#include "tables/St_g2t_vpd_hit_Table.h"
//#include "tables/St_g2t_track_Table.h"
//#include "tables/St_g2t_tpc_hit_Table.h"

typedef vector<StTofMCSlat> tofMCSlatVector;
typedef tofMCSlatVector::iterator tofMCSlatVecIter;


static RanluxEngine engine;

ClassImp(StTofSimMaker)

StTofSimMaker::StTofSimMaker(const char *name):StMaker(name){/* nope */}

StTofSimMaker::~StTofSimMaker(){/* nope */}

Int_t StTofSimMaker::Init(){
  mGeomDb = new StTofGeometry();
  mGeomDb->init();
  mCalibDb = new StTofCalibration();
  mCalibDb->init();
  mSimDb = new StTofSimParam();
  mSimDb->init();
  
  // book histograms
  mdE = new TH1F("energy deposition","dE",100,0.,0.011);
  mdS = new TH1F("distance","ds",100,0.,10);
  mNumberOfPhotoelectrons = new TH1F("number of photoelectrons","nphe",1000,0,5000);
  mT = new TH1F("delay corrected tof","time",100,0.,12e-7);
  mTime = new TH1F("only hit-pos resolution added","tt",100,0.,12e-7);
  mTime1 = new TH1F("fully corrected tof","tt1",100,0.,120e-7);
  mPMlength = new TH1F("distance in slat","length",100,0,22);

  return StMaker::Init();
}


Int_t StTofSimMaker::Make(){
  cout << "StTofSimMaker  Make() starts" << endl;

  StTofSlatCollection *mSlatCollection =  new StTofSlatCollection;
  tofMCSlatVector tofMC;

  // Check for GEANT data
  St_DataSet *geantData = GetInputDS("geant");
  if (geantData) {
    St_DataSetIter geantIter(geantData);

    // TOFp hits
    St_g2t_ctf_hit *g2t_tof_hit = (St_g2t_ctf_hit*) geantIter("g2t_tof_hit");
    if (g2t_tof_hit){
      g2t_ctf_hit_st* tof_hit = g2t_tof_hit->GetTable();
      int numberOfTofHits = g2t_tof_hit->GetNRows();
      cout << "TOF #hits: " << numberOfTofHits << endl;
      
      tofMCSlatVector MCSlatVec;
      MCSlatVec.clear();
      for (int i=0;i<numberOfTofHits;i++,tof_hit++){
	MCSlatVec.push_back(detectorResponse(tof_hit));
      }

      // build ... from MCSlatVec which may have entries with same slatId
      tofMCSlatVector slatTempVec = MCSlatVec;
      tofMCSlatVector slatErasedVec = slatTempVec;
      tofMCSlatVecIter slatTempIter, slatErasedIter;
	
      while (slatTempVec.size()!=0){
	unsigned short fastTdc;
	int nFired=0, accumNPhe=0;
	float accumDe=0., accumDs=0., fastTof=0.;
	slatTempIter=slatTempVec.begin();
	slatErasedIter=slatErasedVec.begin();
	fastTof = slatTempIter->mcInfo().mTof;
	fastTdc = slatTempIter->tdc();
	    
	while(slatErasedIter!=slatErasedVec.end()){
	  if(slatTempIter->slatIndex() == slatErasedIter->slatIndex()){
	    nFired++;
	    accumDe += slatErasedIter->mcInfo().mDe;
	    accumDs += slatErasedIter->mcInfo().mDs;
	    accumNPhe += slatErasedIter->mcInfo().mNPhe;
	    fastTof = min(fastTof, slatErasedIter->mcInfo().mTof);
	    fastTdc = min(fastTdc, slatErasedIter->tdc());
	      
	    slatErasedVec.erase(slatErasedIter);
	    slatErasedIter--;
	  }
	  slatErasedIter++;
	}
	StTofMCSlat MCSlat = *slatTempIter;
	MCSlat.setNHits(nFired);
	MCSlat.setNPhe(accumNPhe);
	MCSlat.setDe(accumDe);
	MCSlat.setDs(accumDs);
	MCSlat.setTof(fastTof);
	MCSlat.setTdc(fastTdc);
	tofMC.push_back(MCSlat);
	
	slatTempVec = slatErasedVec;
      }
      cout << "StTofSimMaker::make()  vector size from " << MCSlatVec.size()
	   << " to " << tofMC.size() << endl;
	
      //////
	
      fillRaw();
      if (mSimDb->elec_noise() < 0)  electronicNoise();
      fillEvent();
      for (unsigned int i=0;i<tofMC.size(); i++){
	StTofMCSlat *MCSlatPtr = new StTofMCSlat();
	*MCSlatPtr = tofMC[i];
	//cout << *MCSlatPtr ;
	mSlatCollection->push_back(MCSlatPtr);
      }
    }
    else
      cout << "StTofSimMaker Make()  no TOF hits found" << endl;

    // pVPD section hits
    St_g2t_vpd_hit *g2t_vpd_hit = (St_g2t_vpd_hit*) geantIter("g2t_vpd_hit");
    if (g2t_vpd_hit){
      //   g2t_vpd_hit_st* vpd_hit = g2t_vpd_hit->GetTable();
      int numberOfVpdHits = g2t_vpd_hit->GetNRows();
      cout << "VPD #hits: " << numberOfVpdHits << endl;
    }
    else
      cout << "StTofSimMaker Make()  no VPD hits found" << endl;    
  }


  // send of to StEvent
  StTofCollection *mTheTofCollection = new StTofCollection();
  for (size_t j=0;j<mSlatCollection->size();j++){
    mTheTofCollection->addSlat(dynamic_cast<StTofMCSlat*>(mSlatCollection->getSlat(j))); 
    // mTheTofCollection->addSlat(mSlatCollection->getSlat(j)); 
}

  mEvent = (StEvent*) GetInputDS("StEvent");
  mEvent->setTofCollection(mTheTofCollection);
  


  // verify existence of tofCollection in StEvent (mEvent) 
  cout << "StTofSimMaker: verifying TOF StEvent data ..." << endl;
  StTofCollection *mmTheTofCollection = mEvent->tofCollection();
  if(mmTheTofCollection) {
    cout << " + StEvent tofCollection Exists" << endl;
    if(mmTheTofCollection->slatsPresent())
      cout << " + StEvent TofSlatCollection Exists" << endl;
    else
      cout << " - StEvent TofSlatCollection DOES NOT Exist" << endl;
    if(mmTheTofCollection->hitsPresent())
      cout << " + StEvent TofHitCollection Exists" << endl;
    else
      cout << " - StEvent TofHitCollection DOES NOT Exist" << endl;
  }
  else {
    cout << " - StEvent tofCollection DOES NOT Exist" << endl;
    cout << " - StEvent TofSlatCollection DOES NOT Exist" << endl;
    cout << " - StEvent TofHitCollection DOES NOT Exist" << endl;
  }

  cout << "StTofSimMaker  Make() finished" << endl;
  return kStOK;
}


StTofMCSlat StTofSimMaker::detectorResponse(g2t_ctf_hit_st* tof_hit)
{


#ifdef TOFSIM_DEBUG
    // dump the g2t structure ...
    cout << setw( 3) << tof_hit->id      << setw( 4) << tof_hit->next_tr_hit_p
	 << setw( 4) << tof_hit->track_p << setw( 8) << tof_hit->volume_id
	 << setw(13) << tof_hit->de      << setw(11) << tof_hit->ds
	 << setw(12) << tof_hit->p[0]    << setw(12) << tof_hit->p[1]
	 << setw(12) << tof_hit->p[2]    << setw( 7) << tof_hit->s_track
	 << setw(13) << tof_hit->tof     << setw(10) << tof_hit->x[0]
	 << setw(10) << tof_hit->x[1]    << setw(10) << tof_hit->x[2]
	 << endl;
#endif
    // skip the consistency checks for now,

    // determine eta and phi indices from hitpoint and slatId
    StThreeVectorD hitPoint  = StThreeVectorD(tof_hit->x[0],
					      tof_hit->x[1],
					      tof_hit->x[2]);
    //fg    StTofCross*  mTofCross = new StTofCross;
    int slatId= (int) mGeomDb->tofSlatCrossId(hitPoint);
    int volId = (int) mGeomDb->tofSlatCrossId(tof_hit->volume_id);

    if (slatId != volId){
      cout << "StTofSimMaker::Make  Warning: volume_id ("<< volId 
           << ") and hit ("<<slatId<<") inconsistent. Switching to volumeid."<< endl;
      slatId=volId;
    }

    // retrieve dbase parameters
    float zmin = mGeomDb->tofSlat(slatId).z_min;
    float zmax = mGeomDb->tofSlat(slatId).z_max;
    float cosang  = mGeomDb->tofSlat(slatId).cosang;
    
    float length = (zmax-tof_hit->x[2])/cosang ;
    float max_distance = (zmax-zmin)/cosang ;
	
    if (length>max_distance || length<0){
      cout <<  "StTofSimMaker:  length="<<length<<" max="<<max_distance
	   << " zmin="<<zmin<<" zmax="<<zmax<<" coasng="<<cosang<<endl;
      mGeomDb->printSlat(slatId);
    }

	
    // do the slat response modelling similar to cts
    long numberOfPhotoelectrons;
    if (mSimDb->slat_para()){
      cout << "StTofSimMaker  Slat Response Table not implemented yet. "
      << " Switching to exponential model instead" <<endl;
    }
    numberOfPhotoelectrons = long (tof_hit->de * slatResponseExp(length));

    // prepare some random generator stuff
    RandGauss random(engine);

    // calculate TOFs with all kinds of resolutions
    float time= tof_hit->tof + length*mSimDb->delay();
    float resl=  mSimDb->time_res() * sqrt(length);
    if (resl<50e-12) resl=50e-12;
    float tt =  tof_hit->tof + (float) random.shoot()* resl;
    float tt1 =  time +  (float) random.shoot()* mSimDb->start_res();

#ifdef TOFSIM_HIST
    // fill the histograms
    mPMlength->Fill(length);
    mdE->Fill(tof_hit->de);
    mdS->Fill(tof_hit->ds);
    mNumberOfPhotoelectrons->Fill(numberOfPhotoelectrons);
    mT->Fill(time);
    mTime->Fill(tt);
    mTime1->Fill(tt1);
#endif

    // fill or update the mcInfo structure of StTofMCslat
    StTofMCSlat slat;
    slat.setSlatIndex(slatId);
    StTofMCInfo slatData;
    slatData.mNHits = 1; // slats will be reorderd and #hits/slat recounted
    slatData.mNPhe  = numberOfPhotoelectrons;
    slatData.mDe    = tof_hit->de;
    slatData.mDs    = tof_hit->ds;
    slatData.mTof =  tof_hit->tof;// warning: this is *NOT* correct ... prop-delay is not included !!!
    slatData.mTime = time;
    slatData.mMTime = tt;
    slatData.mMTimeL = tt1;
    slatData.mPmLength = length;
    slatData.mSLength = tof_hit->s_track;
    StThreeVectorD hitMomentum  = StThreeVectorD(tof_hit->p[0],
						 tof_hit->p[1],
						 tof_hit->p[2]);
    slatData.mPTot = hitMomentum.mag();
    slatData.mGId = tof_hit->id;
    slatData.mTrkId = tof_hit->track_p;
    
    slat.setMCInfo(slatData);

    // do a rough cts-style calibration
    float tdcOffset = mCalibDb->slat(slatId).offset_tdc
                      + random.shoot() * mCalibDb->slat(slatId).ods_tdc;
    unsigned short tdc = (unsigned short)((float)tt/
                              (mCalibDb->slat(slatId).cc_tdc)+ tdcOffset);
    float adcOffset = mCalibDb->slat(slatId).offset_adc
                      + random.shoot() * mCalibDb->slat(slatId).ods_adc;
    unsigned short adc =  (unsigned short)((float)numberOfPhotoelectrons
                                     * mSimDb->nphe_to_adc() + adcOffset);
    slat.setAdc(adc);
    slat.setTdc(tdc);

    cout << "StTofmcInfo updating slatId " << slatId << " " << slatData;

    // this part considers X-talk between slats (based on parameter below) ...
    //    cout << "PHYSNOISE PARAMETER: " << mSimDb->phys_noise() << endl;
  
    return slat;
}



float StTofSimMaker::slatResponseExp(float& dz)
{
  // Exponential model for slat response
  return mSimDb->GeV_2_n_photons() * mSimDb->cath_eff()
         * mSimDb->cath_surf() * mSimDb->surf_loss()
         * exp(-dz/mSimDb->attlen());
}


Int_t StTofSimMaker::Finish(){

#ifdef TOFSIM_HIST
  cout << "StTofSimMaker::Finish  writing tofsim.root ...";
  TFile theFile("tofsim.root","RECREATE","tofsim");
  theFile.cd();
  mdE->Write();
  mdS->Write();
  mNumberOfPhotoelectrons->Write();
  mT->Write();
  mTime->Write();
  mTime1->Write();
  mPMlength->Write();
  cout << "done"<<endl;
#endif
  return kStOK;
}


void StTofSimMaker::fillRaw(){
  //fill the adc and tdc entries.
}

void StTofSimMaker::electronicNoise(){
}


 void StTofSimMaker::fillEvent(){
}
