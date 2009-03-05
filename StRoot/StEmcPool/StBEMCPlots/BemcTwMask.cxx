#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
using namespace std;

#include <TPad.h>
#include <TH2.h>
#include <TF1.h>
#include <TEnv.h>

#include <TPaveStats.h>
#include <TStyle.h> // for gPad
#include <TROOT.h> // for gROOT
#include <TSystem.h>
#include <TMath.h>

#include <StRoot/StEmcPool/StBEMCPlots/BemcTwMask.h>

//--------------------------------------
bool  useBtowMask(const Char_t *fname, BemcTwMask *m) {
  
  cout << "BEMCPlots::useTwMask(\"" << fname << "\") ..." << endl;
 
  string line; 
  Int_t softid, cr, seq, TWmask, HTmask ,TPmask, tp;
  Int_t plot, channelinplot, nok=0;
  Float_t ped;
  Char_t name[1000];

  ifstream BemcStatusFile(fname);

  if (!BemcStatusFile)
    {
      cout << "BemcTwMask:: useBTowMask Failed due to NO FILE"<<endl;
      m->clear();
      return 0;
    }
  
  if (BemcStatusFile)
    {
      if (BemcStatusFile.is_open()) 
	{
	  while (! BemcStatusFile.eof() )//loop until you get to the end of the file
	    {
	      std::getline(BemcStatusFile,line);//get next line in files
	      size_t comment=line.find("#",0); //check for # in line
	      size_t trigger=line.find("trigger",0);//check for "trigger" in line
	      if (comment!=string::npos) continue; //skip comment line
	      if (trigger!=string::npos) continue; //skip trigger mask line
	      istringstream iss; iss.str(line);//parse line for 
	      iss>>name>>softid>>cr>>seq>>TWmask>>HTmask>>TPmask>>ped>>tp;
	      if (softid<1||softid>4800) return 0;//sanity check

	      if (TWmask!=1){//if status !=good
		
		channelinplot=softid;
		
		if ((softid>0)&&(softid<=1220))  plot=0;
		if ((softid>1220)&&(softid<=2400)) plot=1;
		if ((softid>2400)&&(softid<=3540)) plot=2;
		if ((softid>3540)&&(softid<=4800)) plot=3;
		
		int jj=m->crG[plot].GetN();
		m->crG[plot].SetPoint(jj,channelinplot,-10);
		m->crG[plot].SetPoint(jj+1,channelinplot,100000);
		m->crG[plot].SetPoint(jj+2,channelinplot,-10);
		//m->crG2[plot].SetPoint(jj,channelinplot,-10);
		//m->crG2[plot].SetPoint(jj+1,channelinplot,100000);
		//m->crG2[plot].SetPoint(jj+2,channelinplot,-10);
		
		nok++;
	      }
	    }
	}
      BemcStatusFile.close();
    } 
  
  m->nMask=nok;
   
  int i;
  for(i=0;i<m->nPl;i++){
    m->crG[i].SetMarkerStyle(23);
    m->crG[i].SetMarkerColor(kRed);
    //m->crG2[i].SetMarkerStyle(1);
    //m->crG2[i].SetLineColor(kRed);
  }
 
}
