/*!
 * \class StV0Vertex 
 * \author Thomas Ullrich, Sep 1999
 */
/***************************************************************************
 *
 * $Id: StV0Vertex.h,v 2.5 2002/03/08 20:28:37 ullrich Exp $
 *
 * Author: Thomas Ullrich, Sep 1999
 ***************************************************************************
 *
 * Description:
 *
 ***************************************************************************
 *
 * $Log: StV0Vertex.h,v $
 * Revision 2.5  2002/03/08 20:28:37  ullrich
 * Custom Streamer written.
 *
 * Revision 2.4  2002/02/22 22:56:53  jeromel
 * Doxygen basic documentation in all header files. None of this is required
 * for QM production.
 *
 * Revision 2.3  2001/04/05 04:00:46  ullrich
 * Replaced all (U)Long_t by (U)Int_t and all redundant ROOT typedefs.
 *
 * Revision 2.2  2001/03/24 03:35:01  perev
 * clone() -> clone() const
 *
 * Revision 2.1  1999/10/28 22:28:04  ullrich
 * Adapted new StArray version. First version to compile on Linux and Sun.
 *
 * Revision 2.0  1999/10/12 18:43:26  ullrich
 * Completely Revised for New Version
 *
 **************************************************************************/
#ifndef StV0Vertex_hh
#define StV0Vertex_hh
#include "StVertex.h"
#include "StContainers.h"

class dst_v0_vertex_st;

class StV0Vertex : public StVertex {
public:
    StV0Vertex();
    StV0Vertex(const dst_vertex_st&, const dst_v0_vertex_st&);
    // StV0Vertex(const StV0Vertex&);            use default
    // StV0Vertex& operator=(const StV0Vertex&); use default
    ~StV0Vertex();
    
    StVertexId            type() const;
    unsigned int          numberOfDaughters() const;
    StTrack*              daughter(StChargeSign sign);
    const StTrack*        daughter(StChargeSign sign) const;
    StTrack*              daughter(unsigned int);
    const StTrack*        daughter(unsigned int) const;
    StPtrVecTrack         daughters(StTrackFilter&);
    void                  addDaughter(StTrack*);
    void                  removeDaughter(StTrack*);

    float                 dcaDaughterToPrimaryVertex(StChargeSign sign) const;
    float                 dcaDaughters() const;
    float                 dcaParentToPrimaryVertex() const;
    const StThreeVectorF& momentumOfDaughter(StChargeSign sign) const;
    StThreeVectorF        momentum() const;
    
    void setDcaDaughterToPrimaryVertex(StChargeSign, float);
    void setMomentumOfDaughter(StChargeSign, const StThreeVectorF&);
    void setDcaDaughters(float);
    void setDcaParentToPrimaryVertex(float);

private:
    StPtrVecTrack    mDaughters;
    Float_t          mDcaDaughtersToPrimaryVertex[2];
    //    StThreeVectorF   mMomentumOfDaughters[2];
    StThreeVectorF   mMomentumOfDaughters_0; // negative
    StThreeVectorF   mMomentumOfDaughters_1; // positive
    Float_t          mDcaDaughters;
    Float_t          mDcaParentToPrimaryVertex;

    StObject* clone() const;
    ClassDef(StV0Vertex,2)
};
#endif
