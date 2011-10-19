#ifndef L2wEemc2012_H
#define L2wEemc2012_H
/*********************************************************************
 * $Id: L2wEemc2012.h,v 1.2 2011/10/19 15:39:45 jml Exp $
 * \author Jan Balewski, MIT, 2008 
 *********************************************************************
 * Descripion:
 *  example algo finding list of 2x2 BTOW clusters based on ET-seed list produced by L2-btow-calib algo
 * features: uses 2D array (ieta vs. iphi )
 * Limitations:
 *   - ignores seeds at the  edges
 *   - double processing  for neighbours seeds
 *********************************************************************
 */

#ifdef IS_REAL_L2  //in l2-ana  environment
  #include "L2VirtualAlgo2012.h"
#else
  #include "StTriggerUtilities/L2Emulator/L2algoUtil/L2VirtualAlgo2012.h"
#endif

class L2Histo;
class L2EmcGeom;

class L2wEemc2012 : public  L2VirtualAlgo2012 {
 private:

  //..................... params set in initRun
  int   par_dbg;
  float par_EtThresh;

 //............... preserved for Decision(),
  float highestEt[L2eventStream2012::mxToken]; 
  
  // utility methods
  void  createHisto();
  
 public:
  L2wEemc2012(const char* name,const char *uid, L2EmcDb* db, L2EmcGeom *geo, char* outDir, int resOff);
  int   initRunUser( int runNo, int *rc_ints, float *rc_floats);
  void  finishRunUser();// at the end of each run
  void  computeUser(int token);
  bool  decisionUser(int token, int *myL2Result);
};

#endif 

/**********************************************************************
  $Log: L2wEemc2012.h,v $
  Revision 1.2  2011/10/19 15:39:45  jml
  2012

  Revision 1.1  2011/10/18 15:11:45  jml
  adding 2012 algorithms

  Revision 1.1  2009/03/28 19:43:53  balewski
  2009 code



*/

