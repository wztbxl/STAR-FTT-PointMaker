#include "StFwdUtils/StFwdAnalysisMaker.h"
#include "StFwdTrackMaker/Common.h"

#include "TMath.h"

#include <limits>
#include <map>
#include <string>
#include <string>
#include <vector>

#include "StBFChain/StBFChain.h"

#include "StEvent/StEvent.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StHelixModel.h"
#include "StEvent/StPrimaryTrack.h"
#include "StEvent/StRnDHit.h"
#include "StEvent/StRnDHitCollection.h"
#include "StEvent/StTrack.h"
#include "StEvent/StTrackGeometry.h"
#include "StEvent/StTrackNode.h"
#include "StEvent/StPrimaryVertex.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StTrackDetectorInfo.h"
#include "StEvent/StFttPoint.h"
#include "StEvent/StFcsHit.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFttCollection.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StTriggerData.h"
#include "StEvent/StFstHitCollection.h"
#include "StEvent/StFstHit.h"
#include "StEvent/StFwdTrackCollection.h"
#include "StChain/StChainOpt.h"

#include "StEventUtilities/StEventHelper.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrack.h"


#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_vertex_Table.h"
#include "tables/St_g2t_event_Table.h"

#include "StarMagField/StarMagField.h"

#include "St_base/StMessMgr.h"
#include "StarClassLibrary/StPhysicalHelix.hh"
#include "StarClassLibrary/SystemOfUnits.h"


#include "TROOT.h"
#include "TLorentzVector.h"
#include "StEvent/StFwdTrack.h"
#include "StFcsDbMaker/StFcsDb.h"

//________________________________________________________________________
StFwdAnalysisMaker::StFwdAnalysisMaker() : StMaker("fwdAna"){};
int StFwdAnalysisMaker::Finish() { 
    
    auto prevDir = gDirectory;
        
    // output file name
    string name = "FwdAna.root";
    TFile *fOutput = new TFile(name.c_str(), "RECREATE");
    fOutput->cd();
    for (auto nh : mHists) {
        nh.second->SetDirectory(gDirectory);
        nh.second->Write();
    }

    // restore previous directory
    gDirectory = prevDir;

    return kStOk; 
}
//________________________________________________________________________
int StFwdAnalysisMaker::Init() { 
    LOG_DEBUG << "StFwdAnalysisMaker::Init" << endm; 

    mHists["fwdMultFailed"] = new TH1F("fwdMultFailed", ";N_{ch}^{FWD}; counts", 100, 0, 100);
    mHists["fwdMultAll"] = new TH1F("fwdMultAll", ";N_{ch}^{FWD}; counts", 100, 0, 100);
    mHists["fwdMultGood"] = new TH1F("fwdMultGood", ";N_{ch}^{FWD}; counts", 100, 0, 100);
    mHists["fwdMultFST"] = new TH1F("fwdMultFST", ";N_{ch}^{FWD}; counts", 100, 0, 100);
    mHists["fwdMultEcalMatch"] = new TH1F("fwdMultEcalMatch", ";N_{ch}^{FWD}; counts", 100, 0, 100);
    mHists["fwdMultHcalMatch"] = new TH1F("fwdMultHcalMatch", ";N_{ch}^{FWD}; counts", 100, 0, 100);
    mHists["eta"] = new TH1F("eta", ";#eta; counts", 100, 0, 5);
    mHists["phi"] = new TH1F("phi", ";#phi; counts", 100, -3.1415926, 3.1415926);

    mHists["ecalMatchPerTrack"] = new TH1F("ecalMatchPerTrack", ";N_{match} / track; counts", 5, 0, 5);
    mHists["hcalMatchPerTrack"] = new TH1F("hcalMatchPerTrack", ";N_{match} / track; counts", 5, 0, 5);

    mHists["ecalEnergy"] = new TH1F("ecalEnergy", ";Energy; counts", 100, 0, 15);
    mHists["hcalEnergy"] = new TH1F("hcalEnergy", ";Energy; counts", 100, 0, 15);

    mHists["ecalXY"] = new TH2F( "ecalXY", ";ecalX;ecalY", 200, -200, 200, 200, -200, 200 );
    mHists["hcalXY"] = new TH2F( "hcalXY", ";hcalX;hcalY", 200, 0, 50, 200, 0, 50 );

    mHists[ "ecaldX" ] = new TH1F( "ecaldX", ";dx (trk - ecal); counts", 400, -200, 200 );
    mHists[ "matchedEcaldX" ] = new TH1F( "matchedEcaldX", ";dx (trk - ecal); counts", 400, -200, 200 );
    mHists[ "ecaldY" ] = new TH1F( "ecaldY", ";dy (trk - ecal); counts", 400, -200, 200 );
    mHists[ "matchedEcaldY" ] = new TH1F( "matchedEcaldY", ";dy (trk - ecal); counts", 400, -200, 200 );
    mHists[ "ecaldR" ] = new TH1F( "ecaldR", ";dr (trk - ecal); counts", 400, 0, 400 );
    mHists[ "ecalMindR" ] = new TH1F( "ecalMindR", ";dr (trk - ecal); counts", 400, 0, 400 );
    mHists[ "matchedEcaldR" ] = new TH1F( "matchedEcaldR", ";dr (trk - ecal); counts", 400, 0, 400 );

    mHists[ "hcaldX" ] = new TH1F( "hcaldX", ";dx (trk - hcal); counts", 400, -200, 200 );
    mHists[ "matchedHcaldX" ] = new TH1F( "matchedHcaldX", ";dx (trk - hcal); counts", 400, -200, 200 );
    mHists[ "hcaldY" ] = new TH1F( "hcaldY", ";dy (trk - hcal); counts", 400, -200, 200 );
    mHists[ "matchedHcaldY" ] = new TH1F( "matchedHcaldY", ";dy (trk - hcal); counts", 400, -200, 200 );
    mHists[ "hcaldR" ] = new TH1F( "hcaldR", ";dr (trk - hcal); counts", 400, 0, 400 );
    mHists[ "hcalMindR" ] = new TH1F( "hcalMindR", ";dr (trk - hcal); counts", 400, 0, 400 );
    mHists[ "matchedHcaldR" ] = new TH1F( "matchedHcaldR", ";dr (trk - hcal); counts", 400, 0, 400 );

    mHists[ "trkEcalX" ] = new TH2F( "trkEcalX", ";trkX;ecalX", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkEcalY" ] = new TH2F( "trkEcalY", ";trkY;ecalY", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkEcalMinX" ] = new TH2F( "trkEcalMinX", ";trkX;ecalX", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkEcalMinY" ] = new TH2F( "trkEcalMinY", ";trkY;ecalY", 300, -150, 150, 300, -150, 150 );

    mHists[ "trkHcalX" ] = new TH2F( "trkHcalX", ";trkX;hcalX", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkHcalY" ] = new TH2F( "trkHcalY", ";trkY;hcalY", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkHcalMinX" ] = new TH2F( "trkHcalMinX", ";trkX;hcalX", 300, -150, 150, 300, -150, 150 );
    mHists[ "trkHcalMinY" ] = new TH2F( "trkHcalMinY", ";trkY;hcalY", 300, -150, 150, 300, -150, 150 );

    return kStOK;
}
//________________________________________________________________________
int StFwdAnalysisMaker::Make() {
    LOG_DEBUG << "StFwdAnalysisMaker::Make" << endm;
    StEvent *event = (StEvent *)GetDataSet("StEvent");
    if (event){
        StFttCollection *fttCol = event->fttCollection();
        if (fttCol){
            LOG_INFO << "The Ftt Collection has " << fttCol->numberOfPoints() << " points" << endm;
        }
    }
    long long itStart = FwdTrackerUtils::nowNanoSecond();
    ProcessFwdTracks();
    // ProcessFwdMuTracks();
    LOG_DEBUG << "Processing Fwd Tracks took: " << (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e6 << " ms" << endm;
    return kStOK;
} // Make
//________________________________________________________________________
void StFwdAnalysisMaker::Clear(const Option_t *opts) { LOG_DEBUG << "StFwdAnalysisMaker::CLEAR" << endm; }
//________________________________________________________________________
void StFwdAnalysisMaker::ProcessFwdTracks(  ){
    // This is an example of how to process fwd track collection
    LOG_INFO << "StFwdAnalysisMaker::ProcessFwdTracks" << endm;
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    if (!stEvent)
        return;
    StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
    if (!ftc)
        return;

    LOG_INFO << "Checking FcsCollection" << endm;
    StFcsCollection *fcs = stEvent->fcsCollection();
    if (!fcs) return;

    StFcsDb *mFcsDb = static_cast<StFcsDb *>(GetDataSet("fcsDb"));

    // if (ftc->tracks().size() > 4) return;

    size_t fwdMultEcalMatch = 0;
    size_t fwdMultHcalMatch = 0;
    size_t fwdMultFST = 0;

    mHists[ "fwdMultAll" ]->Fill( ftc->tracks().size() );

    size_t nGood = 0;
    size_t nFailed = 0;
    for ( auto fwdTrack : ftc->tracks() ){
        if ( !fwdTrack->didFitConvergeFully() ) {
            nFailed++;
            continue;
        }
        nGood++;
        LOG_INFO << TString::Format("StFwdTrack[ nProjections=%lu, nFTTSeeds=%lu, nFSTSeeds=%lu, mPt=%f ]", fwdTrack->mProjections.size(), fwdTrack->mFTTPoints.size(), fwdTrack->mFSTPoints.size(), fwdTrack->momentum().perp()) << endm;
        
        LOG_INFO << "StFwdTrack has " << fwdTrack->ecalClusters().size() << " ecal matches" << endm;
        LOG_INFO << "StFwdTrack has " << fwdTrack->hcalClusters().size() << " hcal matches" << endm;

        mHists["ecalMatchPerTrack"]->Fill( fwdTrack->ecalClusters().size() );
        mHists["hcalMatchPerTrack"]->Fill( fwdTrack->hcalClusters().size() );
        
        if (fwdTrack->mFSTPoints.size() > 0){
            fwdMultFST ++;
        }

        mHists["eta"]->Fill( fwdTrack->momentum().pseudoRapidity() );
        mHists["phi"]->Fill( fwdTrack->momentum().phi() );
    
        // ecal proj
        float c[9];
        int detId = kFcsWcalId;
        TVector3 ecalXYZ;
        TVector3 ecapP;

        StFwdTrackProjection ecalProj = fwdTrack->getProjectionFor( detId, 0 );
        StFwdTrackProjection hcalProj = fwdTrack->getProjectionFor( kFcsHcalId, 0 );
        LOG_INFO << "EcalProj z= " << ecalProj.mXYZ.z() << endm;
        LOG_INFO << "HcalProj z= " << hcalProj.mXYZ.z() << endm;

        for ( int iEcal = 0; iEcal < fwdTrack->ecalClusters().size(); iEcal++ ){
            StFcsCluster *clu = fwdTrack->ecalClusters()[iEcal];
            LOG_INFO << "Ecal clu detId = " << clu->detectorId() << endm;
            mHists["ecalEnergy"]->Fill( clu->energy() );

            StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(clu->detectorId(), clu->x(), clu->y());
            float dx = ecalProj.mXYZ.x() - xyz.x();
            float dy = ecalProj.mXYZ.y() - xyz.y();
            float dr = sqrt(dx*dx + dy*dy);
            mHists["matchedEcaldX"]->Fill( dx );
            mHists["matchedEcaldY"]->Fill( dy );
            mHists["matchedEcaldR"]->Fill( dr );
        }

        if (ecalProj.mXYZ.z() > 500){
            double mindR = 999;
            StFcsCluster * cclu = nullptr; // closet cluster
            for ( int iDet = 0; iDet < 2; iDet++ ){
                for( int i = 0; i < fcs->clusters(iDet).size(); i++){
                    StFcsCluster * clu = fcs->clusters(iDet)[i];

                    StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(clu->detectorId(), clu->x(), clu->y());
                    mHists["ecalXY"]->Fill( xyz.x(), xyz.y() );

                    float dx = ecalProj.mXYZ.x() - xyz.x();
                    float dy = ecalProj.mXYZ.y() - xyz.y();
                    float dr = sqrt(dx*dx + dy*dy);

                    if ( fabs(dy) < 25 )
                        mHists[ "ecaldX" ]->Fill( dx );
                    if ( fabs(dx) < 25 )
                        mHists[ "ecaldY" ]->Fill( dy );
                    mHists[ "ecaldR" ]->Fill( dr );
                    if ( dr < mindR ){
                        mindR = dr;
                        cclu = clu;
                    }

                    mHists[ "trkEcalX" ] -> Fill( ecalProj.mXYZ.x(), xyz.x() );
                    mHists[ "trkEcalY" ] -> Fill( ecalProj.mXYZ.y(), xyz.y() );

                }
            }
            mHists[ "ecalMindR" ]->Fill( mindR );
            if (cclu){
                StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(cclu->detectorId(), cclu->x(), cclu->y());
                mHists[ "trkEcalMinX" ] -> Fill( ecalProj.mXYZ.x(), xyz.x() );
                mHists[ "trkEcalMinY" ] -> Fill( ecalProj.mXYZ.y(), xyz.y() );
            }
        }

        if (hcalProj.mXYZ.z() > 500){
            
            double mindR = 999;
            StFcsCluster * cclu = nullptr;
            for ( int iDet = 2; iDet < 4; iDet++ ){
                for( int i = 0; i < fcs->clusters(iDet).size(); i++){
                    StFcsCluster * clu = fcs->clusters(iDet)[i];
                    if (!clu) continue;
                    StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(clu->detectorId(), clu->x(), clu->y());
                    mHists["hcalXY"]->Fill( xyz.x(), xyz.y() );

                    float dx = hcalProj.mXYZ.x() - xyz.x();
                    float dy = hcalProj.mXYZ.y() - xyz.y();
                    float dr = sqrt(dx*dx + dy*dy);

                    if ( fabs(dy) < 25 )
                        mHists[ "hcaldX" ]->Fill( dx );
                    if ( fabs(dx) < 25 )
                        mHists[ "hcaldY" ]->Fill( dy );
                    mHists[ "hcaldR" ]->Fill( dr );

                    if ( dr < mindR ){
                        mindR = dr;
                        cclu = clu;
                    }

                    mHists[ "trkHcalX" ] -> Fill( hcalProj.mXYZ.x(), xyz.x() );
                    mHists[ "trkHcalY" ] -> Fill( hcalProj.mXYZ.y(), xyz.y() );
                }
            }
            mHists[ "hcalMindR" ]->Fill( mindR );
            if (cclu){
                StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(cclu->detectorId(), cclu->x(), cclu->y());
                mHists[ "trkHcalMinX" ] -> Fill( hcalProj.mXYZ.x(), xyz.x() );
                mHists[ "trkHcalMinY" ] -> Fill( hcalProj.mXYZ.y(), xyz.y() );
            }
        }

        if (fwdTrack->ecalClusters().size() > 0)
            fwdMultEcalMatch++;
        if (fwdTrack->hcalClusters().size() > 0)
            fwdMultHcalMatch++;

    } // Loop ftc->tracks()

    mHists[ "fwdMultGood" ]->Fill( nGood );
    mHists[ "fwdMultFailed" ]->Fill( nFailed );
    mHists["fwdMultFST"]->Fill( fwdMultFST );
    mHists["fwdMultHcalMatch"]->Fill( fwdMultHcalMatch );
    mHists["fwdMultEcalMatch"]->Fill( fwdMultEcalMatch );
} // ProcessFwdTracks

//________________________________________________________________________
void StFwdAnalysisMaker::ProcessFwdMuTracks(  ){
    // This is an example of how to process fwd track collection
    LOG_INFO << "StFwdAnalysisMaker::ProcessFwdMuTracks" << endm;
    StMuDstMaker *mMuDstMaker = (StMuDstMaker *)GetMaker("MuDst");
    if(!mMuDstMaker) {
        LOG_WARN << " No MuDstMaker ... bye-bye" << endm;
        return;
    }
    StMuDst *mMuDst = mMuDstMaker->muDst();
    if(!mMuDst) {
        LOG_WARN << " No MuDst ... bye-bye" << endm;
        return;
    }
    StMuFwdTrackCollection * ftc = mMuDst->muFwdTrackCollection();
    if (!ftc) return;

    StMuFcsCollection *fcs = mMuDst->muFcsCollection();
    if (!fcs) return;

    cout << "Number of StMuFwdTracks: " << ftc->numberOfFwdTracks() << endl;

    if (ftc->numberOfFwdTracks() > 4 ) return;

    StFcsDb *mFcsDb = static_cast<StFcsDb *>(GetDataSet("fcsDb"));
    

    size_t fwdMultFST = 0;
    size_t fwdMultEcalMatch = 0;
    size_t fwdMultHcalMatch = 0;

    for ( size_t iTrack = 0; iTrack < ftc->numberOfFwdTracks(); iTrack++ ){
        StMuFwdTrack * muFwdTrack = ftc->getFwdTrack( iTrack );
        // LOG_DEBUG << TString::Format("StMuFwdTrack[ nProjections=%lu, nFTTSeeds=%lu, nFSTSeeds=%lu, mPt=%f ]", muFwdTrack->mProjections.size(), muFwdTrack->mFTTPoints.size(), muFwdTrack->mFSTPoints.size(), muFwdTrack->momentum().Pt()) << endm;

        LOG_INFO << "StMuFwdTrack has " << muFwdTrack->mEcalClusters.GetEntries() << " Ecal matched" << endm;
        LOG_INFO << "StMuFwdTrack has " << muFwdTrack->mHcalClusters.GetEntries() << " Hcal matched" << endm;

        mHists["eta"]->Fill( muFwdTrack->momentum().Eta() );
        mHists["phi"]->Fill( muFwdTrack->momentum().Phi() );

        if (muFwdTrack->mFSTPoints.size() > 0){
            fwdMultFST ++;
        }

        if (muFwdTrack->mEcalClusters.GetEntries() > 0)
            fwdMultEcalMatch++;
        if (muFwdTrack->mHcalClusters.GetEntries() > 0)
            fwdMultHcalMatch++;

        
        // ecal proj
        float c[9];
        int detId = kFcsWcalId;
        TVector3 ecalXYZ;
        TVector3 ecapP;

        StMuFwdTrackProjection ecalProj;
        bool foundEcalProj = muFwdTrack->getProjectionFor( detId, ecalProj, 0 );

        if (foundEcalProj){

            // LOG_INFO << "EcalProj z= " << ecalProj.mXYZ.Z() << endm;
            for( int i = 0; i < fcs->numberOfClusters(); i++){
                StMuFcsCluster * clu = fcs->getCluster(i);

                if ( clu->detectorId() > 1 ) continue;

                if ( clu->energy() < 1 ) continue;
                StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(clu->detectorId(), clu->x(), clu->y());

                float dx = ecalProj.mXYZ.X() - xyz.x();
                float dy = ecalProj.mXYZ.Y() - xyz.y();
                float dr = sqrt(dx*dx + dy*dy);

                mHists[ "ecaldX" ]->Fill( dx );
                mHists[ "ecaldY" ]->Fill( dy );
                mHists[ "ecaldR" ]->Fill( dr );

                mHists[ "trkEcalX" ] -> Fill( ecalProj.mXYZ.X(), xyz.x() );

            }

            
        }

        
        for ( size_t i = 0; i < muFwdTrack->mEcalClusters.GetEntries(); i++ ){
            auto c = (StMuFcsCluster*) muFwdTrack->mEcalClusters.At(i);
            if (!c) continue;
            mHists["ecalEnergy"]->Fill( c->energy() );
            
            LOG_INFO << "eCal Cluster detId = " << c->detectorId() << endm;
            StThreeVectorD xyz = mFcsDb->getStarXYZfromColumnRow(c->detectorId(), c->x(), c->y());
            mHists["ecalXY"]->Fill( xyz.x(), xyz.y() );

            if (foundEcalProj)
                mHists["matchedEcaldX"]->Fill( ecalProj.mXYZ.X() - xyz.x() );

            
        }

        mHists["ecalMatchPerTrack"]->Fill( muFwdTrack->mEcalClusters.GetEntries() );
        mHists["hcalMatchPerTrack"]->Fill( muFwdTrack->mHcalClusters.GetEntries() );

        for ( size_t i = 0; i < muFwdTrack->mHcalClusters.GetEntries(); i++ ){
            auto c = (StMuFcsCluster*) muFwdTrack->mHcalClusters.At(i);
            if (!c) continue;
            mHists["hcalEnergy"]->Fill( c->energy() );

            mHists["hcalXY"]->Fill( c->x(), c->y() );
        }

    }
    mHists["fwdMult"]->Fill( ftc->numberOfFwdTracks() );
    mHists["fwdMultFST"]->Fill( fwdMultFST );
    mHists["fwdMultHcalMatch"]->Fill( fwdMultHcalMatch );
    mHists["fwdMultEcalMatch"]->Fill( fwdMultEcalMatch );

}
