/*********************************** 
** Fri Mar 28 12:39:58 1997
** Automatically generated by stafGen 
** DO NOT EDIT THIS FILE 
***********************************/ 
 
/*------------------ INCLUDES ----*/ 
#include <stream.h> 
#include <stdlib.h> 
#include <string.h> 
#include "kuip.h" 
 
/*------------------ ASPS --------*/ 
#include "asuLib.h" 
#include "emlLib.h" 
#include "socLib.h" 
#include "spxLib.h" 
#include "tdmLib.h" 
#include "levLib.h" 
#include "duiLib.h" 
#include "dioLib.h" 
#include "amiLib.h" 
#include "tbrLib.h" 
 
/*------------------ PAMS --------*/ 
 
/*------------------ PROTOTYPES --*/ 
extern CC_P int ami_load(amiBroker *broker); 
extern CC_P void unknown(); 
extern CC_P int stafArgs(int argc, char **argv); 
#ifdef TNT 
extern CC_P void staf_paw_init_();  
#else 
extern CC_P void staf_kuip_init_();  
#endif 
extern CC_P void staf_banner(FILE* stream); 
 
/*================== MAIN ========*/ 
int main(int argc, char** argv) 
{ 
   stafArgs(argc,argv); 
#ifdef TNT 
   staf_paw_init_();  
#else 
   staf_kuip_init_();  
#endif 
 
/*------------------ START -------*/ 
#ifdef ASU 
   asu_init(); asu_start(); 
#endif /* ASU */ 
 
#ifdef EML 
   eml_init(); eml_start(); 
#endif /* EML */ 
 
#ifdef SOC 
   soc_init(); soc_start(); 
#endif /* SOC */ 
 
#ifdef SPX 
   spx_init(); spx_start(); 
#endif /* SPX */ 
 
#ifdef TDM 
   tdm_init(); tdm_start(); 
#endif /* TDM */ 
 
#ifdef DUI 
   dui_init(); dui_start(); 
#endif /* DUI */ 
 
#ifdef LEV 
   lev_init(); lev_start(); 
#endif /* LEV */ 
 
#ifdef DIO 
   dio_init(); dio_start(); 
#endif /* DIO */ 
 
#ifdef AMI 
   ami_init(); ami_start(); 
#endif /* AMI */ 
 
#ifdef TBR 
   tbr_init(); tbr_start(); 
#endif /* TBR */ 
 
 
/*------------------ INPUT -------*/ 
   staf_banner(stdout); 
   printf("STAF - ready for user input.\n"); 
   ku_what(unknown); 
   printf("STAF - done with user input.\n"); 
 
/*------------------ STOP --------*/ 
#ifdef TBR 
   tbr_stop(); 
#endif /* TBR */ 
 
#ifdef AMI 
   ami_stop(); 
#endif /* AMI */ 
 
#ifdef DIO 
   dio_stop(); 
#endif /* DIO */ 
 
#ifdef DUI 
   dui_stop(); 
#endif /* DUI */ 
 
#ifdef LEV 
   lev_stop(); 
#endif /* LEV */ 
 
#ifdef TDM 
   tdm_stop(); 
#endif /* TDM */ 
 
#ifdef SPX 
   spx_stop(); 
#endif /* SPX */ 
 
#ifdef SOC 
   soc_stop(); 
#endif /* SOC */ 
 
#ifdef EML 
   eml_stop(); 
#endif /* EML */ 
 
#ifdef ASU 
   asu_stop(); 
#endif /* ASU */ 
 
 
   return 1; 
} 
 
/*================== UNKNOWN =====*/ 
void unknown(){} 
 
/*================== AMI_LOAD ====*/ 
int ami_load(amiBroker *broker) 
{ 
 
   return TRUE; 
} 
 
