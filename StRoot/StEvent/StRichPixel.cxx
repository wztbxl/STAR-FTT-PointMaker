/***************************************************************************
 *
 * $Id: StRichPixel.cxx,v 2.3 2000/01/10 17:12:21 lasiuk Exp $
 *
 * Author: Thomas Ullrich, July 1999
 ***************************************************************************
 *
 * Description:
 *
 ***************************************************************************
 *
 * $Log: StRichPixel.cxx,v $
 * Revision 2.3  2000/01/10 17:12:21  lasiuk
 * remove dst_rch_pixel dependency;
 * change stored data to a single long;
 * modify unpacking routines;
 *
 * Revision 2.2  1999/10/28 22:26:22  ullrich
 * Adapted new StArray version. First version to compile on Linux and Sun.
 *
 * Revision 2.1  1999/10/13 19:45:05  ullrich
 * Initial Revision
 *
 **************************************************************************/
#include "StRichPixel.h"

static const char rcsid[] = "$Id: StRichPixel.cxx,v 2.3 2000/01/10 17:12:21 lasiuk Exp $";

ClassImp(StRichPixel)

StRichPixel::StRichPixel()
    : mPackedData(0)
{ /* nopt */ }

StRichPixel::StRichPixel(ULong_t rawData)
   : mPackedData(rawData)
{ /* noop */ }

StRichPixel::~StRichPixel() { /* noop */ }
    
Int_t
StRichPixel::operator==(const StRichPixel& p) const
{
    return (p.mPackedData  == mPackedData);
}

Int_t
StRichPixel::operator!=(const StRichPixel& p) const
{
    return !(*this == p);  // use operator==()
}

UShort_t
StRichPixel::pad() const
{
    return (mPackedData & 0xff);  // first 8 bits
}

UShort_t
StRichPixel::row() const
{
    return ( (mPackedData>>8) & 0xff);  // second 8 bits
}

UShort_t
StRichPixel::adc()  const
{
    return ( (mPackedData>>16) & 0x3ff); // 10bits
}
