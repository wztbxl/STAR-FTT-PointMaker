/***************************************************************************
 *
 * $Id: StSvtEmbeddingMaker.h,v 1.7 2004/07/01 13:54:29 caines Exp $
 *
 * Author: Selemon Bekele
 ***************************************************************************
 *
 * Description: Svt Embedding Maker class
 *
 ***************************************************************************
 *
 * $Log: StSvtEmbeddingMaker.h,v $
 * Revision 1.7  2004/07/01 13:54:29  caines
 * Changes to the simulation maker from the review
 *
 * Revision 1.5  2004/02/24 15:53:21  caines
 * Read all params from database
 *
 * Revision 1.4  2004/01/22 16:30:47  caines
 * Getting closer to a final simulation
 *
 * Revision 1.3  2003/11/30 20:51:48  caines
 * New version of embedding maker and make OnlSeqAdj a stand alone maker
 *
 * Revision 1.2  2003/09/10 19:47:37  perev
 * ansi corrs
 *
 * Revision 1.1  2003/07/31 19:18:09  caines
 * Petrs improved simulation code
 *
 **************************************************************************/
#ifndef STAR_StSvtEmbeddingMaker
#define STAR_StSvtEmbeddingMaker

#ifndef StMaker_H
#include "StMaker.h"
#endif


class StSvtHybridData;
class StSvtData;
class StSvtHybridCollection;
class StSvtHybridPixelsD;
/*!
 *
 * \class  StSvtEmbeddingMaker
 * \author Chaloupka
 * \date   2004/07/29
 * \brief  Automaticaly takes care of creation of the background in SVT Slow Simulator.
 *         This is embedding or only simulated noise.
 *
 */
class StSvtEmbeddingMaker : public StMaker
{

public:
  
  StSvtEmbeddingMaker(const char* name = "SvtEmbedding");
  virtual ~StSvtEmbeddingMaker();
  
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();
  virtual Int_t InitRun(int runumber);
    
  void setBackGround(Bool_t backgr,double backgSigma);      //default is TRUE and 1.8
  void setDoEmbedding(Bool_t doIt);                       //this allows(when set to FALSE) to force the embedding maker to ignore raw data 
  void SetPedRmsPreferences(Bool_t usePixelRMS, Bool_t useHybridRMS);                    // allows to disable reading RMS from database
							  //and create simple background -just a plain simulation; default is TRUE
  
private:
  void ReadPedRMSfromDb();
  void GetSvtData();
  void GetPedRMS();
  void ClearMask();
  void AddRawData();
  void CreateBackground();
  double  MakeGaussDev(double sigma);

  StSvtData*               mSimPixelColl;   //!
  StSvtData*               mRealDataColl;       //!
  StSvtHybridCollection*   mPedColl;            //!
  StSvtHybridCollection*   mPedRMSColl;         //!

  ///parameterss - read from database 
  double mBackGSigma;  ///default value if individiual RMS are not available 
    Bool_t mBackGrOption; //shell I generate background
  Bool_t mDoEmbedding;  ///shell I try to run embedding if posssible 
  Bool_t mUsePixelRMS;  ///try to read individual pixel RMS from database, default TRUE
  Bool_t mUseHybridRMS; ///try to read individual hybrid RMS from database, default TRUE
  
  Bool_t mRunningEmbedding;   /// can I realy run embedding - ie. missing DAQ maker?
  Bool_t mMask[128*240];

  //global variables for temporary store in the loop
  StSvtHybridPixelsD  *mCurrentPixelData;
  int mCurrentIndex;
  
  ClassDef(StSvtEmbeddingMaker,1)

};

#endif


