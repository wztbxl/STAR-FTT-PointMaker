/*:>--------------------------------------------------------------------
**: FILE:       track_propagator.c.template
**: HISTORY:
**:             00jan93-v000a-hpl- Created by stic Version
**:  Id: idl.y,v 1.14 1998/08/28 21:22:28 fisyak Exp  
**:<------------------------------------------------------------------*/
#include "track_propagator.h"

long track_propagator_(
  TABLE_HEAD_ST         *gtrack_h,      DST_TRACK_ST           *gtrack ,
  TABLE_HEAD_ST         *target_h,  EGR_PROPAGATE_ST           *target )
{
/*:>--------------------------------------------------------------------
**: ROUTINE:    track_propagator_
**: DESCRIPTION: Physics Analysis Module ANSI C template.
**:             This is an ANSI C Physics Analysis Module template
**:             automatically generated by stic from track_propagator.idl.
**:             Please edit comments and code.
**: AUTHOR:     hpl - H.P. Lovecraft, hplovecraft@cthulhu.void
**: ARGUMENTS:
**:       IN:
**:             target    - PLEASE FILL IN DESCRIPTION HERE
**:            target_h   - header Structure for target
**:    INOUT:
**:             gtrack    - PLEASE FILL IN DESCRIPTION HERE
**:            gtrack_h   - header Structure for gtrack
**:      OUT:
**: RETURNS:    STAF Condition Value
**:>------------------------------------------------------------------*/
  
  
  long  i, iflag, iiflag=1;
  float bfld=5., GRADDEG=57.2958;
  float psi, pt, tanl, x0, y0, z0, xp[2], xout[4], p[3], xv[3], xx0[3];
  float trk[6], r1, r2, xc1[2], xc2[2], x[2], y[2], cut, x1[3];
  long  q;

  for(i=0; i<gtrack_h->nok; i++) 
    {
      psi  = gtrack[i].psi;
      q    = (long) gtrack[i].icharge;
      tanl = gtrack[i].tanl;
      if(gtrack[i].invpt != 0 )
	{
	  pt    = 1/gtrack[i].invpt ;
	}
      else
	{
	  /*	    message('Error, invpt); */
	  pt = 0.01;
	}
      x0   = gtrack[i].x0 ;
      y0   = gtrack[i].y0 ;
      z0   = gtrack[i].z0;
      trk[0] = gtrack[i].x0;
      trk[1] = gtrack[i].y0;
      trk[2] = gtrack[i].z0;
      trk[3] = gtrack[i].psi;
      trk[4] = gtrack[i].tanl;
      trk[5] = gtrack[i].icharge;
      r1     = ((float) gtrack[i].icharge)*(1/gtrack[i].invpt)/(0.0003*bfld);
     
      
      if(target->iflag == 1 || target->iflag == 2)
	{
	  xp[0]= target->x[0]; 
	  xp[1]= target->x[1];
	  
	  prop_project_track_(&psi, &q, &pt, &tanl, &x0, &y0, &bfld, &z0, xp, xout); 
	  
	  gtrack[i].x0 = xout[0];
	  gtrack[i].y0 = xout[1];
	  gtrack[i].z0 = xout[2];
	  gtrack[i].psi= xout[3];

	  if (target->iflag == 2)
	    {
	      xv[0] = target->x[0];
	      xv[1] = target->x[1];
	      xv[2] = target->x[2];
	      x1[0] = xout[0];
	      x1[1] = xout[1];
	      x1[2] = xout[2];
	      trk[0] = gtrack[i].x0;
	      trk[1] = gtrack[i].y0;
	      trk[2] = gtrack[i].z0;
	      trk[3] = gtrack[i].psi;
	      prop_track_mom_(trk, &r1, &iiflag, p, &bfld);  
	      prop_fine_approach_(xv, x1, p, &xx0); 
	      gtrack[i].x0 = xx0[0];
	      gtrack[i].y0 = xx0[1];
	      gtrack[i].z0 = xx0[2];
	    }
	  return STAFCV_OK;
	}
      
      if(target->iflag == 3)
	{
	  prop_circle_param_(trk, xc1, &r1, &bfld); 
	  xc2[0] = xc2[1] = 0.;
	  r2 = target->r;
	  cut = 0.4;
	  prop_vzero_geom_(&cut, xc1, xc2, &r1, &r2, x, y, &iflag); 
	  if(iflag == 5)
	    {
	      /* message("Can't project to the circle target->r! "); */
	      continue;
	    }
	  else if (iflag == 3)
	    {
      	      xp[0] = x[0];
	      xp[1] = y[0];
	    }
	  else
	    {
	      if(  ((x0-x[0])*(x0-x[0]) + (y0-y[0])*(y0-y[0])) <
		   ((x0-x[1])*(x0-x[1]) + (y0-y[1])*(y0-y[1]))   )
		{ 
		  xp[0] = x[0];
		  xp[1] = y[0];
		}
	      else
		{
		  xp[0] = x[1];
		  xp[1] = y[1];
		} 
	    }
	  
	  prop_project_track_(&psi, &q, &pt, &tanl, &x0, &y0, &bfld, &z0, xp, xout); 
	  
	  gtrack[i].x0 = xout[0];
	  gtrack[i].y0 = xout[1];
	  gtrack[i].z0 = xout[2];
	  gtrack[i].psi= xout[3];
	}
    }
}
  













