#ifndef StiTrackNode_H
#define StiTrackNode_H 1

#include <iostream.h>
#include <stdlib.h>

#include "StiDefaultMutableTreeNode.h"
#include "Messenger.h"

class StiHit;
class StiDetector;

class StiTrackNode : public StiDefaultMutableTreeNode
{
public:
    
    void reset();
    void set(int depth, StiHit * h);
    void setAsCopyOf(const StiTrackNode * node);
    void     setHit(StiHit * h)  {  
	hit = h;  
	//    if(hit!=NULL){ detector = NULL; }
	if(hit!=0){ detector = 0; }
    }
    StiHit * getHit() const      {  return hit;}
    friend ostream& operator<<(ostream& os, const StiTrackNode& n);
    
    const StiDetector *getDetector() const;
    void setDetector(const StiDetector *pDetector);
    
    void setDedx(double e) {dedx=e;}
    double getDedx() const {return dedx;}
    
protected:   
    
    //StiTrackNode(): hit(NULL), detector(NULL){}
    StiTrackNode(): hit(0), detector(0){
      if(s_pMessenger==NULL){ 
        s_pMessenger = Messenger::instance(MessageType::kNodeMessage);
      }
    }
    
    double        dedx;
    StiHit      * hit;  
    const StiDetector * detector; // used if not hit for node

    static Messenger *s_pMessenger;
};

#endif



