/***************************************************************************
 *
 * $Id: StHbtParticleCut.h,v 1.3 2000/01/07 23:21:17 laue Exp $
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: StHbtMaker package
 *     base class for particle-wise cuts
 * Note:    Users DO NOT inherit from this class!
 *          The base classes StHbtTrackCut and StHbtV0Cut inherit from this,
 *          and users inherit from those
 * 
 ***************************************************************************
 *
 * $Log: StHbtParticleCut.h,v $
 * Revision 1.3  2000/01/07 23:21:17  laue
 * 0.) all 'ClassDef(...)' put inside #ifdef __ROOT__  #endif
 * 1.) unnecessary includes of 'StMaker.h' deleted
 *
 * Revision 1.2  1999/12/03 22:24:34  lisa
 * (1) make Cuts and CorrFctns point back to parent Analysis (as well as other way). (2) Accommodate new PidTraits mechanism
 *
 * Revision 1.1  1999/10/15 01:56:50  lisa
 * Important enhancement of StHbtMaker - implement Franks CutMonitors
 * ----------------------------------------------------------
 * This means 3 new files in Infrastructure area (CutMonitor),
 * several specific CutMonitor classes in the Cut area
 * and a new base class in the Base area (StHbtCutMonitor).
 * This means also changing all Cut Base class header files from .hh to .h
 * so we have access to CutMonitor methods from Cint command line.
 * This last means
 * 1) files which include these header files are slightly modified
 * 2) a side benefit: the TrackCuts and V0Cuts no longer need
 * a SetMass() implementation in each Cut class, which was stupid.
 * Also:
 * -----
 * Include Franks StHbtAssociationReader
 * ** None of these changes should affect any user **
 *
 * Revision 1.3  1999/09/17 22:37:59  lisa
 * first full integration of V0s into StHbt framework
 *
 * Revision 1.2  1999/07/06 22:33:19  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.1.1.1  1999/06/29 16:02:56  lisa
 * Installation of StHbtMaker
 *
 **************************************************************************/


#ifndef StHbtParticleCut_hh
#define StHbtParticleCut_hh

#include "StHbtMaker/Infrastructure/StHbtTypes.hh"
#include "StHbtMaker/Infrastructure/StHbtCutMonitorHandler.h"

class StHbtParticleCut : public StHbtCutMonitorHandler {

public:

  StHbtParticleCut(){/* no-op */};   // default constructor. - Users should write their own
  virtual ~StHbtParticleCut(){/* no-op */};  // destructor

  virtual StHbtString Report() =0;    // user-written method to return string describing cuts

  double Mass(){return mMass;};       // mass of the particle being selected
  virtual void SetMass(const double& mass) {mMass = mass;};

  virtual StHbtParticleType Type()=0;


  // the following allows "back-pointing" from the CorrFctn to the "parent" Analysis
  friend class StHbtAnalysis;
  StHbtAnalysis* HbtAnalysis(){return myAnalysis;};


protected:
  double mMass;
  //  StHbtParticleType mType;            // tells whether cut is on Tracks or V0's
  StHbtAnalysis* myAnalysis;

#ifdef __ROOT__
  ClassDef(StHbtParticleCut, 0)
#endif
};

#endif
