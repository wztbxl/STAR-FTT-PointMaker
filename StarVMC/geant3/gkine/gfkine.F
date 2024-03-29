*
* $Id: gfkine.F,v 1.1.1.2 2009/02/01 17:10:30 fisyak Exp $
*
* $Log: gfkine.F,v $
* Revision 1.1.1.2  2009/02/01 17:10:30  fisyak
* Merge with .DEV2
*
* Revision 1.1.1.1  2008/12/10 20:48:28  fisyak
* Merge with macos version
*
* Revision 1.1.1.1  2005/05/25 22:37:53  fisyak
* Alice version of geant3 (-minicern)
*
* Revision 1.2  2003/11/28 11:23:56  brun
* New version of geant321 with all geant3 routines renamed from G to G3
*
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:40  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.2  1998/02/10 16:04:18  japost
*   Comments are correctly redistributed.
*
* Revision 1.1.1.1  1995/10/24 10:21:18  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE G3FKINE(ITRA,VERT,PVERT,IPART,NVERT,UBUF,NWBUF)
*
************************************************************************
*                                                                      *
*           Retrieves Vertex and Track parameters                      *
*           -------------------------------------                      *
*                                                                      *
*  Retrieves parameters for a track with a long lifetime.              *
*                                                                      *
*  Input:                                                              *
*    ITRA      track number for which parameters are requested         *
*                                                                      *
*  Output:                                                             *
*    VERT      vector origin of the track                              *
*    PVERT     4 momentum components at the track origin               *
*    IPART     particle type (=0 if track ITRA does not exist)         *
*    NVERT     vertex number origin of the track                       *
*    UBUF      user words stored in GSKINE.                            *
*                                                                      *
*    ==>Called by : <USER>, GTREVE                                     *
*       Author    R.Brun  *********                                    *
*                                                                      *
************************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gctrak.inc"
      DIMENSION VERT(3),PVERT(4),UBUF(*)
C.
C.    ------------------------------------------------------------------
C.
      IPART  = 0
      IF (JKINE.LE.0)GO TO 99
      IF (ITRA.LE.0)GO TO 99
      IF (ITRA.GT.NTRACK)GO TO 99
C
      JK = LQ(JKINE- ITRA)
      IF (JK.LE.0)GO TO 99
      NV = Q(JK + 6)
      IF (JVERTX.LE.0)GO TO 99
      JV = LQ(JVERTX- NV)
      IF (JV.LE.0)GO TO 99
      DO 10 I=1,3
  10  VERT(I) = Q(JV + I)
*     TOFG=Q(JV+4)
C
      DO 20 I=1,4
  20  PVERT(I) = Q(JK + I)
      IPART = Q(JK + 5)
      NVERT = Q(JK + 6)
      NLINK = IQ(JK-2)
      NWBUF  = 0
      IF (NLINK.LE.0)GO TO 99
      JKU = LQ(JK- 1)
      IF (JKU.LE.0)GO TO 99
      NWBUF = IQ(JKU-1)
      DO 30 I=1,NWBUF
  30  UBUF(I)=Q(JKU+I)
C
  99  RETURN
      END
