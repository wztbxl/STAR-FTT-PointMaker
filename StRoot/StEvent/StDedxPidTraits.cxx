/***************************************************************************
 *
 * $Id: StDedxPidTraits.cxx,v 2.3 1999/11/16 14:11:38 ullrich Exp $
 *
 * Author: Thomas Ullrich, Sep 1999
 ***************************************************************************
 *
 * Description:
 *
 ***************************************************************************
 *
 * $Log: StDedxPidTraits.cxx,v $
 * Revision 2.3  1999/11/16 14:11:38  ullrich
 * Changed variance to sigma.
 *
 * Revision 2.5  1999/11/29 17:07:24  ullrich
 * Moved method() from StTrackPidTraits to StDedxPidTraits.cxx
 *
 * Revision 2.4  1999/11/23 15:56:23  ullrich
 * Added clone() method. Was pure virtual.
 *
 * Revision 2.3  1999/11/16 14:11:38  ullrich
 * Changed variance to sigma.
 *
 * Revision 2.2  1999/10/28 22:25:01  ullrich
 * Adapted new StArray version. First version to compile on Linux and Sun.
 *
 * Revision 2.1  1999/10/13 19:44:31  ullrich
 * Initial Revision
 *
 **************************************************************************/
#include "tables/St_dst_dedx_Table.h"
#include "StDedxPidTraits.h"

ClassImp(StDedxPidTraits)

    mNumberOfPoints(0), mDedx(0), mSigma(0) { /* noop */ }

StDedxPidTraits::StDedxPidTraits() :
    mNumberOfPoints(0), mDedx(0), mSigma(0), mMethod(0) { /* noop */ }
    StTrackPidTraits(det, meth),
    mNumberOfPoints(n), mDedx(dedx), mSigma(sig) { /* noop */ }
                                 UShort_t n, Float_t dedx, Float_t sig) :
    StTrackPidTraits(det),
    mNumberOfPoints(n), mDedx(dedx), mSigma(sig), mMethod(meth) { /* noop */ }
    mNumberOfPoints(t.ndedx), mDedx(t.dedx[0]), mSigma(t.dedx[1]) { /* noop */ }
    StTrackPidTraits(t),
    mNumberOfPoints(t.ndedx), mDedx(t.dedx[0]),
    mSigma(t.dedx[1]), mMethod(t.method){ /* noop */ }

StDedxPidTraits::~StDedxPidTraits() { /* noop */ }

UShort_t
StDedxPidTraits::numberOfPoints() const { return mNumberOfPoints; }

Float_t
StDedxPidTraits::sigma() const { return mSigma; }
    }
}

