#ifndef STSVTONLINESEQADJSIM_HH
#define STSVTONLINESEQADJSIM_HH

#include "Stiostream.h"
#include <stdlib.h>
#include "StMaker.h"

class StSvtConfig;
class StSvtHybridPixelsC;
class StSvtHybridPixelsD;
//class StSvtHybridData;
class StSvtHybridCollection;
class StSvtData;
class StSvtDaq;

//used bad anode list is the same as in StSvtSeqAdjMaker => no difference between online and offline


class StSvtOnlineSeqAdjSimMaker:public StMaker
{
public:
  StSvtOnlineSeqAdjSimMaker(const char* name = "SvtOnlineSeqAdj");
  ~StSvtOnlineSeqAdjSimMaker();

  
  void SetKillBadAnodes(bool doit){mKillBadAnodes=doit;}
  bool GetKillBadAnodes(){return mKillBadAnodes;}
  void SetSaveAnode2Raw(bool doit){mSaveAnode2Raw=doit;}
  bool GetSaveAnode2Raw(){return mSaveAnode2Raw;}
  void SetSaveAnode239Raw(bool doit){mSaveAnode239Raw=doit;}
  bool GetSaveAnode239Raw(){return mSaveAnode239Raw;}
  
  void SetNumberTBinsToClear(int n){mNumberTBinsToClear=n;}
  int GetNumberTBinsToClear(){return mNumberTBinsToClear;}
  void SetExtraPixelsBefore(int num){mExtraBefore=num;}
  int GetExtraPixelsBefore(){return mExtraBefore;}
  void SetExtraPixelsAfter(int num){mExtraAfter=num;}
  int GetExtraPixelsAfter(){return mExtraAfter;}

  void Set_n_seq_lo(int num){m_n_seq_lo=num;}
  int Get_n_seq_lo(){return m_n_seq_lo;}
  void Set_n_seq_hi(int num){m_n_seq_hi=num;}
  int Get_n_seq_hi(){return m_n_seq_hi;}
  void Set_thresh_lo(int num){m_thresh_lo=num;}
  int Get_thresh_lo(){return m_thresh_lo;}
  void Set_thresh_hi(int num){m_thresh_hi=num;}
  int Get_thresh_hi(){return m_thresh_hi;}
  void SetPedOffset(int num){mPedOffset=num;}
  int GetPedOffset(){return mPedOffset;}
  
  void SetAdjParams(int thresh_lo,int n_seq_lo,int thresh_hi,int n_seq_hi);
 
 //inherited maker routines
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void  Clear(const char *opt);
  virtual Int_t Finish();
  virtual Int_t InitRun(int runumber); //caled when run number changes

private:
  StSvtData                  *mPixelColl;       //! the simulated data - created for each run InitRun{in beginAnalyses} 
  StSvtData                  *m8bitPixelColl;   //! simulated final result written to 8 bits
  StSvtData                  *mRawData;         //!
  StSvtHybridCollection      *mSvtBadAnodes;    //!
  StSvtConfig                *mConfig;          //!
  StSvtDaq                   *mDaq;             //!

  bool                       mKillBadAnodes;    //!
  int                        mNumberTBinsToClear;//!number of first time bins which are set to 0 by DAQ - default 2
  bool                       mSaveAnode2Raw;     //!
  bool                       mSaveAnode239Raw;   //!

  //number of extra pixels around the sequence to save
  int   mExtraBefore;//!
  int   mExtraAfter;//!

  //parameters for seq adjusting
  int m_n_seq_lo;//!
  int m_n_seq_hi;//!
  int m_thresh_lo;//!
  int m_thresh_hi;//!
  int mPedOffset;//!
  Bool_t mMask[128*240];//!

  //global variables for temporary store in the loop
  StSvtHybridPixelsD  *mCurrentPixelData;//!
  StSvtHybridPixelsC  *mCurrent8bitPixelData;//!
  int mCurrentIndex;

  Int_t GetConfig();
  Int_t GetDaqParams();
  void  GetBadAnodes();
  void  GetPixelData();
  void  SetRawData();

  void  Conversion10to8bit();
  void  KillBadAnodes();
  void  ClearMask();
  void  RawAnodes();
  void  ClearFirstTbins();
  void  SequenceSearch();
  void  WriteMask();
  void  WriteSequence(int anode,int begins, int ends, int NumOfHigh);
  void  FillRawData();
  
 ClassDef(StSvtOnlineSeqAdjSimMaker,1)
};

#endif
