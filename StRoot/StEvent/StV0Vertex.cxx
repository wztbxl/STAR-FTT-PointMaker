/***************************************************************************
 *
 * $Id: StV0Vertex.cxx,v 2.6 2002/03/08 22:16:20 jeromel Exp $
 *
 * Author: Thomas Ullrich, Sep 1999
 ***************************************************************************
 *
 * Description:
 *
 ***************************************************************************
 *
 * $Log: StV0Vertex.cxx,v $
 * Revision 2.6  2002/03/08 22:16:20  jeromel
 * Extra version needed to be read (Thomas)
 *
 * Revision 2.5  2002/03/08 20:28:36  ullrich
 * Custom Streamer written.
 *
 * Revision 2.4  2001/04/05 04:00:59  ullrich
 * Replaced all (U)Long_t by (U)Int_t and all redundant ROOT typedefs.
 *
 * Revision 2.3  2001/03/24 03:35:00  perev
 * clone() -> clone() const
 *
 * Revision 2.2  1999/12/21 15:09:21  ullrich
 * Modified to cope with new compiler version on Sun (CC5.0).
 *
 * Revision 2.1  1999/10/28 22:28:01  ullrich
 * Adapted new StArray version. First version to compile on Linux and Sun.
 *
 * Revision 2.0  1999/10/12 18:43:23  ullrich
 * Completely Revised for New Version
 *
 **************************************************************************/
#include <algorithm>
#include "StV0Vertex.h"
#include "StTrack.h"
#include "StTrackGeometry.h"
#include "tables/St_dst_vertex_Table.h"
#include "tables/St_dst_v0_vertex_Table.h"
#if !defined(ST_NO_NAMESPACES)
using std::fill_n;
using std::copy;
#endif

ClassImp(StV0Vertex)

static const char rcsid[] = "$Id: StV0Vertex.cxx,v 2.6 2002/03/08 22:16:20 jeromel Exp $";

StV0Vertex::StV0Vertex()
{
    mType = kV0VtxId;
    mDaughters.resize(2);
    mDaughters[negative] = 0;
    mDaughters[positive] = 0;
    fill_n(mDcaDaughtersToPrimaryVertex, 2, 0);
    mDcaDaughters             = 0;
    mDcaParentToPrimaryVertex = 0;
}

StV0Vertex::StV0Vertex(const dst_vertex_st& vtx, const dst_v0_vertex_st& v0vtx) : StVertex(vtx)
{
    mType = kV0VtxId;
    mDaughters.resize(2);
    mDaughters[negative] = 0;
    mDaughters[positive] = 0;
    mDcaDaughtersToPrimaryVertex[negative] = v0vtx.dcan;
    mDcaDaughtersToPrimaryVertex[positive] = v0vtx.dcap;
    mMomentumOfDaughters_0.setX(v0vtx.neg_px);
    mMomentumOfDaughters_0.setY(v0vtx.neg_py);
    mMomentumOfDaughters_0.setZ(v0vtx.neg_pz);
    mMomentumOfDaughters_1.setX(v0vtx.pos_px);
    mMomentumOfDaughters_1.setY(v0vtx.pos_py);
    mMomentumOfDaughters_1.setZ(v0vtx.pos_pz);
    mDcaDaughters = v0vtx.dcapn;
    mDcaParentToPrimaryVertex = v0vtx.dcav0;
}

StV0Vertex::~StV0Vertex() { /* noop */ }

StObject*
StV0Vertex::clone() const { return new StV0Vertex(*this); }

StVertexId
StV0Vertex::type() const { return kV0VtxId; }

unsigned int
StV0Vertex::numberOfDaughters() const { return 2; }

StTrack*
StV0Vertex::daughter(unsigned int i)
{
    return i < 2 ?  mDaughters[i] : 0;
}

const StTrack*
StV0Vertex::daughter(unsigned int i) const
{
    return i < 2 ?  mDaughters[i] : 0;
}

StPtrVecTrack
StV0Vertex::daughters(StTrackFilter& filter)
{
    StPtrVecTrack vec;
    for (int i=0; i<2; i++)
        if (filter(mDaughters[i])) vec.push_back(mDaughters[i]);
    return vec;
}

void
StV0Vertex::addDaughter(StTrack* track)
{
    if (track) {
        if (track->geometry()->charge() > 0)
            mDaughters[positive] = track;
        else
            mDaughters[negative] = track;
    }
}

void
StV0Vertex::removeDaughter(StTrack* t)
{
    if (t == mDaughters[positive]) mDaughters[positive] = 0;
    if (t == mDaughters[negative]) mDaughters[negative] = 0;
}

StTrack*
StV0Vertex::daughter(StChargeSign sign)
{
    return mDaughters[sign];
}

const StTrack*
StV0Vertex::daughter(StChargeSign sign) const
{
    return mDaughters[sign];
}

float
StV0Vertex::dcaDaughterToPrimaryVertex(StChargeSign sign) const
{
    return mDcaDaughtersToPrimaryVertex[sign];
}

const StThreeVectorF&
StV0Vertex::momentumOfDaughter(StChargeSign sign) const
{
    return (sign == negative ? mMomentumOfDaughters_0 : mMomentumOfDaughters_1);
}

StThreeVectorF
StV0Vertex::momentum() const
{
    return (mMomentumOfDaughters_0 +
            mMomentumOfDaughters_1);
}

float
StV0Vertex::dcaDaughters() const { return mDcaDaughters; }

float
StV0Vertex::dcaParentToPrimaryVertex() const { return mDcaParentToPrimaryVertex; }

void
StV0Vertex::setDcaDaughterToPrimaryVertex(StChargeSign sign, float val)
{
    mDcaDaughtersToPrimaryVertex[sign] = val;
}

void
StV0Vertex::setMomentumOfDaughter(StChargeSign sign, const StThreeVectorF& v)
{
    if (sign == negative)
	mMomentumOfDaughters_0 = v;
    else
	mMomentumOfDaughters_1 = v;
}

void
StV0Vertex::setDcaDaughters(float val) { mDcaDaughters = val; }

void
StV0Vertex::setDcaParentToPrimaryVertex(float val) { mDcaParentToPrimaryVertex = val; }


/// Custom streamer
void StV0Vertex::Streamer(TBuffer &R__b)
{
    // Stream an object of class .

    if (R__b.IsReading()) {
	UInt_t R__s, R__c;
	Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
	if (R__v > 1) {
	    Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
	    return;
	}

	// *** HACK ***
	//==== process old versions which were saved properly 
	//     by sheer luck.
	StVertex::Streamer(R__b);
	mDaughters.Streamer(R__b);
	R__b >> mDcaDaughtersToPrimaryVertex[0];
	R__b >> mDcaDaughtersToPrimaryVertex[1];
	
	// Old format : had an extra version saved.
	// Was previously saved as StThreeVectorF mMomentumOfDaughters[2];
	// with StThreeVectorF having a custom streamer (therefore the
	// extra version).
	UInt_t R__s2, R__c2;
	Version_t R__v2 = R__b.ReadVersion(&R__s2, &R__c2);
	mMomentumOfDaughters_0.Streamer(R__b);
	mMomentumOfDaughters_1.Streamer(R__b);

	
	// Come back to other data members
	R__b >> mDcaDaughters;
	R__b >> mDcaParentToPrimaryVertex;
	R__b.CheckByteCount(R__s, R__c, Class());
	//====end of old versions
    }
    else {
	Class()->WriteBuffer(R__b,this);
    }
} 



    

    
