//StiStEventFiller.h
/***************************************************************************
 *
 * $Id: StiStEventFiller.h,v 1.2 2002/03/28 04:29:49 calderon Exp $
 *
 * Author: Manuel Calderon de la Barca Sanchez, Mar 2002
 ***************************************************************************
 *
 * $Log: StiStEventFiller.h,v $
 * Revision 1.2  2002/03/28 04:29:49  calderon
 * First test version of Filler
 * Currently fills only global tracks with the following characteristics
 * -Flag is set to 101, as most current global tracks are.  This is not strictly correct, as
 *  this flag is supposed to mean a tpc only track, so really need to check if the track has
 *  svt hits and then set it to the appropriate flag (501 or 601).
 * -Encoded method is set with bits 15 and 1 (starting from bit 0).  Bit 1 means Kalman fit.
 *  Bit 15 is an as-yet unused track-finding bit, which Thomas said ITTF could grab.
 * -Impact Parameter calculation is done using StHelix and the primary vertex from StEvent
 * -length is set using getTrackLength, which might still need tweaking
 * -possible points is currently set from getMaxPointCount which returns the total, and it is not
 *  what we need for StEvent, so this needs to be modified
 * -inner geometry (using the innermostHitNode -> Ben's transformer -> StPhysicalHelix -> StHelixModel)
 * -outer geometry, needs inside-out pass to obtain good parameters at outermostHitNode
 * -fit traits, still missing the probability of chi2
 * -topology map, filled from StuFixTopoMap once StDetectorInfo is properly set
 *
 * This version prints out lots of messages for debugging, should be more quiet
 * when we make progress.
 *
 **************************************************************************/
#ifndef StiStEventFiller_HH
#define StiStEventFiller_HH

//Doxygen class header...
/*! \class StiStEventFiller
    StiStEventFiller is a utilitity class meant to properly convert StiTrack
    objects into StTrack (Global/Primary) objects and hang these on the StEvent
    Track-node.

    \author Manuel Calderon de la Barca Sanchez (Yale Software)
    \note 
    \warning This is a development version, proceed at your own risk!
 */

class StEvent;
class StiTrackContainer;

class StiStEventFiller
{
public:
    StiStEventFiller();
    virtual ~StiStEventFiller();

    ///Fill the event from the track store.
    StEvent* fillEvent(StEvent*, StiTrackContainer*);
    //void operator()(const KalmanTrackMapValType);
    void fillDetectorInfo(StTrackDetectorInfo* detInfo, const StiTrack* kTrack);
    void fillGeometry(StTrack* track, const StiTrack* kTrack, bool outer);
    //void fillTopologyMap(StTrack* track, const StiTrack* kTrack);
    void fillFitTraits(StTrack* track, const StiTrack* kTrack);
    void fillTrack(StTrack* track, const StiTrack* kTrack);
    unsigned short encodedStEventFitPoints(const StiTrack* kTrack); 
    float impactParameter(const StiTrack* kTrack);
private:
    StEvent* mEvent;
    StiTrackContainer* mTrackStore;
};

#endif
