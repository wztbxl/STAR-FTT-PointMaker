// $Id: StBFChain.h,v 1.16 2000/11/25 23:20:48 fisyak Exp $
// $Log: StBFChain.h,v $
// Revision 1.16  2000/11/25 23:20:48  fisyak
// Add dEdx maker to production chain
//
// Revision 1.15  2000/08/06 19:10:15  fisyak
// Clean up BFC, add Bfc table to runco with chain parameters
//
// Revision 1.14  2000/07/26 14:09:03  fisyak
// Split and move emc, remove SetDataBases
//
// Revision 1.13  2000/06/16 23:34:16  fisyak
// replace arrays in parser by TObjArrays
//
// Revision 1.12  2000/05/20 01:03:49  perev
// Sequential event IventNumber added
//
// Revision 1.11  2000/03/14 01:12:52  fisyak
// Split chain files
//
// Revision 1.10  2000/03/09 16:18:35  fisyak
// Make chain table local to the chain
//
// Revision 1.9  1999/11/29 21:38:29  fisyak
// Add Dave Hardtke corrections, clean up print outs
//
// Revision 1.8  1999/11/07 02:26:22  fisyak
// Clean ups
//
// Revision 1.7  1999/11/04 22:21:25  fisyak
// Reorganize chain as Table
//
// Revision 1.6  1999/10/28 01:57:17  fisyak
// Add ParseString
//
// Revision 1.5  1999/10/12 23:13:31  fisyak
// Add AddBefore and AddAfter methods
//
// Revision 1.4  1999/09/24 01:22:51  fisyak
// Reduced Include Path
//
// Revision 1.3  1999/09/12 23:02:43  fisyak
// Add closing xdf and TFile
//
// Revision 1.2  1999/09/08 00:14:06  fisyak
// Add kReverseField option
//
// Revision 1.1  1999/09/02 16:14:43  fine
// new StBFChain library has been introduced to break dependences
//
// Revision 1.4  1999/08/31 00:26:42  fisyak
// Add TFile to BFChain
//
// Revision 1.3  1999/08/10 17:10:52  fisyak
// Exprot EChainOptions into rootcint
//
// Revision 1.2  1999/08/06 14:26:38  fisyak
// put back xdf out option
//
// Revision 1.1  1999/07/29 01:05:23  fisyak
// move bfc to StBFChain
//
#ifndef STAR_StBFChain
#define STAR_StBFChain

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// StBFChain                                                            //
//                                                                      //
// Class to control "BFC" chain                                         //
//                                                                      //
// This class :                                                         //
//   - Initialises the run default parameters                           //
//   - Provides API to Set/Get run parameters                           //
//   - Creates the support lists (TClonesArrays) for the Event structure//
//   - Creates the physics objects makers                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "StChain.h"
#include "TFile.h"
#include "tables/St_Bfc_Table.h"
//_____________________________________________________________________

class St_XDFFile;
class StFileI;
class TObjArray;
#if 0
struct BfcItem {
  Char_t       *Key;      // nick name
  Char_t       *Name;     // maker name
  Char_t       *Chain;    // its chain
  Char_t       *Opts;     // required options
  Char_t       *Maker;    // required Makers
  Char_t       *Libs;     // libraries to be loaded
  Char_t       *Comment;  
  Bool_t        Flag;     // F/T to use it in chain
};
#endif
class StBFChain : public StChain {
 private:
  St_XDFFile          *fXdfOut;  //! xdf output file if any
  TFile               *fTFile;   //TFile associated with the chain
  Bfc_st              *fBFC;      // Private chain
  StFileI             *fSetFiles; //
  TString             *fInFile;   //
  TString             *fFileOut;  //
  TString             *fXdfFile;  //
  
 public:
                       StBFChain(const char *name="bfc");
   virtual            ~StBFChain();
   virtual Int_t       Make(int number){ SetIventNumber(number); return StMaker::Make(number);};
   virtual Int_t       Load();      // *MENU*
   virtual Int_t       Instantiate();      // *MENU*
   virtual Int_t       AddAB (const Char_t *after="",const StMaker *maker=0,const Int_t Opt=1); 
   virtual Int_t       AddAfter  (const Char_t *after, const StMaker *maker) {return AddAB (after,maker);} 
   virtual Int_t       AddBefore (const Char_t *before,const StMaker *maker) {return AddAB (before,maker,-1);} 
   virtual Int_t       ParseString (const TString &tChain, TObjArray &Opt);
   void                SetFlags(const Char_t *Chain="gstar tfs"); // *MENU*
   void                Set_IO_Files(const Char_t *infile=0, const Char_t *outfile=0); // *MENU
   void                SetInputFile(const Char_t *infile=0);                          // *MENU
   void                SetGC(const Char_t *queue=
			     "-s;dst runc;"                     // list of components needed
			     "-q;numberOfPrimaryTracks>1500;"   // example of user query
			     );                                                       // *MENU

   void                SetOutputFile(const Char_t *outfile=0);                        // *MENU
   virtual Int_t       kOpt(const TString *Tag) const; 
   virtual Int_t       kOpt(const Char_t  *Tag) const; 
   virtual void        SetXdfOut(St_XDFFile *xdf=0) {fXdfOut = xdf;}
   virtual void        SetDbOptions();
   virtual void        SetGeantOptions();
   virtual void        SetTreeOptions();
   virtual void        SetOption(const Int_t k);
   virtual void        SetOption(const Char_t*  Opt) {SetOption(kOpt(Opt));}
   virtual void        SetOption(const TString* Opt) {SetOption(kOpt(Opt));}
   virtual void        SetOptionOff(const Char_t*  Opt) {SetOption(-kOpt(Opt));}
   virtual void        SetOptionOff(const TString* Opt) {SetOption(-kOpt(Opt));}
   virtual void        SetTFile(TFile *m) {fTFile = m;}
   virtual Int_t       Finish();
   virtual TFile      *GetTFile() {return fTFile;}
   virtual St_XDFFile *GetXdfOut() {return fXdfOut;}
   virtual Bool_t      GetOption(const Int_t k);
   virtual Bool_t      GetOption(const TString *Opt) {return GetOption(kOpt(Opt));}
   virtual Bool_t      GetOption(const Char_t *Opt)  {return GetOption(kOpt(Opt));}
   virtual const char *GetCVS() const
  {static const char cvs[]="Tag $Name:  $ $Id: StBFChain.h,v 1.16 2000/11/25 23:20:48 fisyak Exp $ built "__DATE__" "__TIME__ ; return cvs;}
   ClassDef(StBFChain, 0)   //StBFChain control class
};
#endif
