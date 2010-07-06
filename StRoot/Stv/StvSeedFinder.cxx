#include "StvSeedFinder.h"
#include "TSystem.h"
#include "StvDraw.h"
#include "StvHit.h"
#include "vector"
#include "StarVMC/GeoTestMaker/StTGeoHelper.h"

ClassImp(StvSeedFinder)
StvSeedFinder* StvSeedFinder::fgSeedFinder=0;
//_____________________________________________________________________________
 StvSeedFinder::StvSeedFinder(const char *name):TNamed(name,"")
 { fDraw=0;fDoShow=0;fgSeedFinder=this;}
//_____________________________________________________________________________
void StvSeedFinder::Clear(const char*)
{
 fSeedHits.clear();
 if(fDraw) fDraw->Clear();
}
//_____________________________________________________________________________
void StvSeedFinder::Show()
{
  if (!fDraw) return;
  fDraw->Trak(fHelix,fSeedHits,kGlobalTrack);
  fDraw->UpdateModified();
}
//_____________________________________________________________________________
void StvSeedFinder::ShowRest()
{
   if (!fDraw) fDraw = NewDraw();
   std::vector<float> myHits;  
   const StVoidArr *hitArr =  StTGeoHelper::Inst()->GetSeedHits();
   int nHits =  hitArr->size();
   for (int iHit=0;iHit<nHits;iHit++) {
     StvHit *stiHit = (StvHit*)(*hitArr)[iHit];
     if (stiHit->timesUsed()) continue;
     const float *x = stiHit->x_g();
    fDraw->Point(x[0],x[1],x[2],kUnusedHit);
   }

   fDraw->UpdateModified();
}
//_____________________________________________________________________________
StvDraw *StvSeedFinder::NewDraw()
{
   StvDraw *dr = new StvDraw();
   dr->SetBkColor(kWhite);
// dr->Style(kUsedHit).Siz()  *=2;
   dr->Style(kUnusedHit).Siz()*=20;
   return dr;
}
//_____________________________________________________________________________
void StvSeedFinder::DoShow(int lev)
{
  fDoShow = lev;
  if (fDoShow) {if (!fDraw) fDraw=NewDraw();}
  else         { delete fDraw;fDraw=0      ;}

}



