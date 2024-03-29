*
* $Id: gfhead.F,v 1.1.1.2 2009/02/01 17:10:18 fisyak Exp $
*
* $Log: gfhead.F,v $
* Revision 1.1.1.2  2009/02/01 17:10:18  fisyak
* Merge with .DEV2
*
* Revision 1.1.1.1  2008/12/10 20:46:57  fisyak
* Merge with macos version
*
* Revision 1.1.1.1  2005/05/25 22:37:25  fisyak
* Alice version of geant3 (-minicern)
*
* Revision 1.2  2003/11/28 11:23:55  brun
* New version of geant321 with all geant3 routines renamed from G to G3
*
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:09  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE G3FHEAD (NWHEAD,IHEAD,NWBUF,UBUF)
C.
C.    ******************************************************************
C.    *                                                                *
C     *       Retrieves the parameter of the Header bank               *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author     M.Maire                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
*
      DIMENSION IHEAD(1),UBUF(1)
*
      IF (JHEAD.NE.0) THEN
         NWHEAD = IQ(JHEAD-1)
         DO 10 I=1,NWHEAD
            IHEAD(I) = IQ(JHEAD+I)
   10    CONTINUE
*
         NWBUF = 0
         IF (LQ(JHEAD-1).NE.0) THEN
            JHU = LQ(JHEAD-1)
            NWBUF = IQ(JHU-1)
            DO 20 I=1,NWBUF
               UBUF(I) = Q(JHU+I)
   20       CONTINUE
         ENDIF
*
      ELSE
         NWHEAD = 0
         WRITE(CHMAIL,1000) IEVENT
         CALL GMAIL(0,0)
      ENDIF
*
 1000 FORMAT(' **** GFHEAD: No Header bank in memory for ievent:',I7)
*
      END
