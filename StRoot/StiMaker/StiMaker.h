//StiMaker.h

#ifndef StiMaker_HH
#define StiMaker_HH

#include <string>
using std::string;

#include "StMaker.h"
#include "StEvent/StEnumerations.h"
#include "Sti/StiKalmanTrackFinder.h"
#include "Sti/StiKalmanTrackFitter.h"

class StEvent;
class StiHit;
class StiTrack;
class StiStEventFiller;
class StiTrackContainer;
class StiTrackFinder;
class StiKalmanTrackFinder;
class StiKalmanTrackFitter;
class StiKalmanTrackNode;
class StiKalmanTrack;
class StiTrackMerger;
class StiToolkit;
class StiMakerParameters;
class StiVertexFinder;
class EventDisplay;
class StiResidualCalculator;
template<class FILTERED> class EditableFilter;


class StiMaker : public StMaker 
{
 public:
    
    StiMaker(const char* name = "StiMaker");
    virtual ~StiMaker();
    virtual void  Clear(const char* opt="");
    virtual Int_t Init();
    virtual Int_t InitDetectors();
    virtual Int_t InitRun(int);
    virtual Int_t Make();
    virtual Int_t Finish();

    virtual const char* GetCVS() const
    {static const char cvs[]="Tag $Name:  $ $Id: StiMaker.h,v 2.20 2006/03/09 22:45:49 didenko Exp $ built "__DATE__" "__TIME__; return cvs;}	

    void setParameters(StiMakerParameters * pars);
    StiMakerParameters * getParameters();

    void setEventDisplay(EventDisplay* eventDisplay) 
      { _eventDisplay = eventDisplay;}
    EventDisplay* getEventDisplay() 
      { return _eventDisplay; }
    StiToolkit * getToolkit();
		void load(TDataSet *);

private:

		double runField;
    StiMakerParameters * _pars;
    bool                 eventIsFinished;
    bool                 _initialized;
    StiToolkit  *        _toolkit;
    StiHitLoader<StEvent,StiDetectorBuilder> * _hitLoader;
    StiTrackFinder *       _seedFinder;
    StiKalmanTrackFinder * _tracker;
    StiKalmanTrackFitter * _fitter;
    StiStEventFiller *    _eventFiller;
    StiTrackContainer *   _trackContainer;
    StiVertexFinder*      _vertexFinder;
    StiResidualCalculator * _residualCalculator;
    EventDisplay *        _eventDisplay;
    EditableFilter<StiTrack> * _loaderTrackFilter;
    EditableFilter<StiHit>   * _loaderHitFilter;
    
    ClassDef(StiMaker,0)
};

//inlines

inline StiToolkit * StiMaker::getToolkit()
{
  return _toolkit;
}
#endif
