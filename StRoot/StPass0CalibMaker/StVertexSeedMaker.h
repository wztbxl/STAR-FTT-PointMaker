/*!
 * \class StVertexSeedMaker 
 * \author G. Van Buren, BNL
 * \version $Id: StVertexSeedMaker.h,v 1.2 2002/01/27 01:56:14 genevb Exp $
 *
 * calculates mean primary vertex positions from
 * suitable events to use as seeds in finding better       
 * primary vertex positions (helpful for low               
 * multiplicity events like pp collisions)                 
 */


#ifndef STAR_StVertexSeedMaker
#define STAR_StVertexSeedMaker

#ifndef StMaker_H
#include "StMaker.h"
#endif
class TNtuple;
class St_vertexSeed;


class StVertexSeedMaker : public StMaker {
 public: 
                  StVertexSeedMaker(const char *name="TpcT0");
   virtual       ~StVertexSeedMaker();
   virtual Int_t Init();
   virtual Int_t Make();
   virtual void PrintInfo();
   virtual void Clear(Option_t *option);
   virtual Int_t Finish();
   virtual void fitData();
   int GetValidityDate();
   int GetValidityTime();
   void UseEventDateTime();
   void UseFillDateTime();
   St_vertexSeed* vertexSeedTable();
   void WriteTableToFile();     //Write drift velocity table (assumes correct trigger offset)
   void SetMinEntries(int entries);  //minimum number of valid events for seed
   void SetMaxX0Err(float err);  //maximum allowed error for x0 
   void SetMaxY0Err(float err);  //maximum allowed error for y0 
   void WriteHistFile();       // Write out vertexseedhist.root file w/results
   void HistFileByDefault();   // Write out file on Finish
   void SetVertexZmax(float zmax);  //Set max z vertex for seed calculation
   void SetVertexZmin(float zmin);  //Set min z vertex for seed calculation
   void SetVertexR2max(float r2max);  //Set max r^2 vertex for seed calculation
   virtual const char *GetCVS() const {
     static const char cvs[]="Tag $Name:  $ $Id: StVertexSeedMaker.h,v 1.2 2002/01/27 01:56:14 genevb Exp $ built "__DATE__" "__TIME__ ;
     return cvs;
   }

 protected:
   Int_t FillAssumed();
   void FillDateTime();

 private:
  TH1F* xdist;
  TH1F* ydist;
  TH1F* xerr;
  TH1F* yerr;
  TNtuple* resNtuple;
  float xguess;
  float yguess;
  float x0_assumed;
  float dxdz_assumed;
  float y0_assumed;
  float dydz_assumed;
  float zvertex;
  float yvertex;
  float xvertex;
  float mult;
  float eventNumber;
  float HIST_MIN;
  float HIST_MAX;
  float zVertexMax; //maximum allowed z vertex for mean calculation
  float zVertexMin; //minimum allowed z vertex for mean calculation
  float r2VertexMax; //minimum allowed radius^2 vertex for mean calculation
  int    fill;
  int    date;
  int    time;
  int    minEntries;
  float    maxX0Err;
  float    maxY0Err;
  Bool_t   mHistOut;
  Bool_t   useEventDateTime;
  double p[4];
  double ep[4];
  double chi;
  double weight;

  ClassDef(StVertexSeedMaker,0)
};


inline void StVertexSeedMaker::UseEventDateTime() {useEventDateTime = kTRUE;}
inline void StVertexSeedMaker::UseFillDateTime() {useEventDateTime = kFALSE;}
inline void StVertexSeedMaker::SetMinEntries(int entries){minEntries = entries; }
inline void StVertexSeedMaker::SetMaxX0Err(float err){maxX0Err = err;}
inline void StVertexSeedMaker::SetMaxY0Err(float err){maxY0Err = err;}
inline int  StVertexSeedMaker::GetValidityDate(){return date;}
inline int  StVertexSeedMaker::GetValidityTime(){return time;}
inline void StVertexSeedMaker::HistFileByDefault(){mHistOut = kTRUE;} 
inline void StVertexSeedMaker::SetVertexZmax(float zmax){zVertexMax = zmax;}
inline void StVertexSeedMaker::SetVertexZmin(float zmin){zVertexMin = zmin;}
inline void StVertexSeedMaker::SetVertexR2max(float r2max){r2VertexMax = r2max;}

#endif
