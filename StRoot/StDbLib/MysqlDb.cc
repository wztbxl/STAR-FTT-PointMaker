/***************************************************************************
 *
 * $Id: MysqlDb.cc,v 1.28 2003/09/18 00:51:59 porter Exp $
 *
 * Author: Laurent Conin
 ***************************************************************************
 *
 * Description: Mysql - SQL Query handler
 *
 ***************************************************************************
 *
 * $Log: MysqlDb.cc,v $
 * Revision 1.28  2003/09/18 00:51:59  porter
 * initialized mhasBinaryQuery flag
 *
 * Revision 1.27  2003/09/16 22:44:17  porter
 * got rid of all ostrstream objects; replaced with ostringstream+string.
 * modified rules.make and added file stdb_streams.h for standalone compilation
 *
 * Revision 1.26  2003/09/02 17:57:49  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.25  2003/07/02 18:39:23  porter
 * added server version to connection notification and host:port to connection failure notification
 *
 * Revision 1.24  2003/04/11 22:47:27  porter
 * Added a fast multi-row write model specifically needed by the daqEventTag
 * writer. Speed increased from about 100Hz to ~3000Hz.  It is only invoked if
 * the table is marked as Non-Indexed (daqTags & scalers). For non-indexed tables
 * which include binary stored data (we don't have any yet), the fast writer  has
 * to invoke a slower buffer so that the rates are a bit slower (~500Hz at 50 rows/insert).
 *
 * Revision 1.23  2003/02/12 22:12:45  porter
 * moved warning message about null columns (checked in 2 days ago) from the
 * depths of the mysql coding into the StDbTable code. This suppresses confusing
 * warnings from tables that have had elements removed but their storage columns
 * still exist in the database.
 *
 * Revision 1.22  2003/02/11 03:22:07  porter
 * added warning message individual columns return null data
 *
 * Revision 1.21  2002/04/16 19:44:49  porter
 * changed non-dbname in arguement of mysql_real_connect from 0 to NULL.  mysql
 * perfers this.  Updated rules.make for local mysql installation on linux
 *
 * Revision 1.20  2002/03/22 19:05:38  porter
 * #-of-retries on server connect increased to 7 with timeout period doubled per
 * retry starting at 1 sec.  DOES NOT work (is ignored) on STAR's Redhat 6.2
 * version of mysqlclient but does on Redhat 7.2. Needed for maintainable
 * multiple mirror servers using dns for round-robin load balancing.
 *
 * Revision 1.19  2002/01/30 15:40:47  porter
 * changed limits on flavor tag & made defaults retrieving more readable
 *
 * Revision 1.18  2001/12/21 04:54:33  porter
 * sped up table definition for emc and changed some ostrstream usage for
 * insure tests
 *
 * Revision 1.17  2001/04/25 17:13:19  perev
 * HPcorrs
 *
 * Revision 1.16  2001/03/31 15:03:46  porter
 * fix bug in StDbManagerImpl::getDbName accidently introduced yesterday
 * & added new diagnostic message in MysqlDb
 *
 * Revision 1.15  2001/03/30 18:48:26  porter
 * modified code to keep Insure from wigging-out on ostrstream functions.
 * moved some messaging into a StDbSql method.
 *
 * Revision 1.14  2001/03/22 19:39:17  porter
 * make a check to avoid mysql bug on linux where selecting an
 * unknown database hangs the connection
 *
 * Revision 1.13  2001/02/09 23:06:24  porter
 * replaced ostrstream into a buffer with ostrstream creating the
 * buffer. The former somehow clashed on Solaris with CC5 iostream (current .dev)
 *
 * Revision 1.12  2001/01/22 18:37:49  porter
 * Update of code needed in next year running. This update has little
 * effect on the interface (only 1 method has been changed in the interface).
 * Code also preserves backwards compatibility so that old versions of
 * StDbLib can read new table structures.
 *  -Important features:
 *    a. more efficient low-level table structure (see StDbSql.cc)
 *    b. more flexible indexing for new systems (see StDbElememtIndex.cc)
 *    c. environment variable override KEYS for each database
 *    d. StMessage support & clock-time logging diagnostics
 *  -Cosmetic features
 *    e. hid stl behind interfaces (see new *Impl.* files) to again allow rootcint access
 *    f. removed codes that have been obsolete for awhile (e.g. db factories)
 *       & renamed some classes for clarity (e.g. tableQuery became StDataBaseI
 *       and mysqlAccessor became StDbSql)
 *
 * Revision 1.11  2000/08/15 22:51:51  porter
 * Added Root2DB class from Masashi Kaneta
 * + made code more robust against requesting data from non-existent databases
 *
 * Revision 1.10  2000/07/27 01:59:18  porter
 * fixed bug in delete vs delete []
 * fixed up LDFLAGS for linux/g++
 *
 * Revision 1.9  2000/06/02 13:37:36  porter
 * built up list of minor changes:
 *  - made buffer more robust for certain null inputs
 *  - fixed small leak in StDbTables & restructure call to createMemory
 *  - added dbRhic as a database domain in StDbDefs
 *  - added setUser() in StDbManager
 *  - added more diagnostic printouts in mysqlAccessor.cc
 *
 * Revision 1.8  2000/03/28 17:03:18  porter
 * Several upgrades:
 * 1. configuration by timestamp for Conditions
 * 2. query by whereClause made more systematic
 * 3. conflict between db-stored comments & number lists resolved
 * 4. ensure endtime is correct for certain query falures
 * 5. dbstl.h->handles ObjectSpace & RogueWave difference (Online vs Offline)
 *
 * Revision 1.7  2000/03/01 20:56:15  porter
 * 3 items:
 *    1. activated reConnect for server timeouts
 *    2. activated connection sharing; better resource utilization but poorer
 *       logging
 *    3. made rollback method in mysqlAccessor more robust (affects writes only)
 *
 * Revision 1.6  2000/02/18 16:58:08  porter
 * optimization of table-query, + whereClause gets timeStamp if indexed
 *  + fix to write multiple rows algorithm
 *
 * Revision 1.5  2000/02/15 20:27:43  porter
 * Some updates to writing to the database(s) via an ensemble (should
 * not affect read methods & haven't in my tests.
 *  - closeAllConnections(node) & closeConnection(table) method to mgr.
 *  - 'NullEntry' version to write, with setStoreMode in table;
 *  -  updated both StDbTable's & StDbTableDescriptor's copy-constructor
 *
 * Revision 1.4  2000/01/27 05:54:31  porter
 * Updated for compiling on CC5 + HPUX-aCC + KCC (when flags are reset)
 * Fixed reConnect()+transaction model mismatch
 * added some in-code comments
 *
 * Revision 1.3  1999/12/07 21:25:25  porter
 * some fixes for linux warnings
 *
 * Revision 1.2  1999/09/30 02:05:59  porter
 * add StDbTime to better handle timestamps, modify SQL content (mysqlAccessor)
 * allow multiple rows (StDbTable), & Added the comment sections at top of
 * each header and src file
 *
 **************************************************************************/
#include "MysqlDb.h"
#include "StDbManager.hh" // for now & only for getting the message service
#include "stdb_streams.h"

//#include "errmsg.h"

#ifdef HPUX
#define freeze(i) str()
#endif


#define CR_MIN_ERROR            2000    /* For easier client code */
#define CR_MAX_ERROR            2999
#define CR_UNKNOWN_ERROR        2000
#define CR_SOCKET_CREATE_ERROR  2001
#define CR_CONNECTION_ERROR     2002
#define CR_CONN_HOST_ERROR      2003
#define CR_IPSOCK_ERROR         2004
#define CR_UNKNOWN_HOST         2005
#define CR_SERVER_GONE_ERROR    2006
#define CR_VERSION_ERROR        2007
#define CR_OUT_OF_MEMORY        2008
#define CR_WRONG_HOST_INFO      2009
#define CR_LOCALHOST_CONNECTION 2010
#define CR_TCP_CONNECTION       2011
#define CR_SERVER_HANDSHAKE_ERR 2012
#define CR_SERVER_LOST          2013
#define CR_COMMANDS_OUT_OF_SYNC 2014
#define CR_NAMEDPIPE_CONNECTION 2015
#define CR_NAMEDPIPEWAIT_ERROR 2016
#define CR_NAMEDPIPEOPEN_ERROR 2017
#define CR_NAMEDPIPESETSTATE_ERROR 2018

#define __CLASS__ "MysqlDb"


static const char* binaryMessage = {"Cannot Print Query with Binary data"};


////////////////////////////////////////////////////////////////////////

MysqlDb::MysqlDb(): mdbhost(0), mdbName(NULL), mdbuser(0), mdbpw(0), mdbPort(0),mdbServerVersion(0),mlogTime(false) {

mhasConnected=false;
mhasBinaryQuery=false;
mtimeout=1;
mQuery=0;
mQueryLast=0;
mRes= new MysqlResult;
 for(int i=0;i<200;i++)cnames[i]=0;

}
//////////////////////////////////////////////////////////////////////

MysqlDb::~MysqlDb(){
if(mQuery) delete [] mQuery;
if(mQueryLast) delete [] mQueryLast;
Release();
if(mRes) delete mRes;
if(mhasConnected)mysql_close(&mData);
if(mdbhost) delete [] mdbhost;
if(mdbuser) delete [] mdbuser;
if(mdbpw)   delete [] mdbpw;
if(mdbName)  delete [] mdbName;
 if(mdbServerVersion) delete [] mdbServerVersion;

}

//////////////////////////////////////////////////////////////////////// 
bool MysqlDb::reConnect(){
#define __METHOD__ "reConnect()"

  bool connected=false;
  unsigned int timeOutConnect=mtimeout;
  while(!connected && timeOutConnect<600){ 
    mysql_options(&mData,MYSQL_OPT_CONNECT_TIMEOUT,(const char*)&timeOutConnect);
    if(mysql_real_connect(&mData,mdbhost,mdbuser,mdbpw,mdbName,mdbPort,NULL,0))
    connected=true;
    if(!connected){
      timeOutConnect*=2;
      ostringstream wm;
      wm<<" Connection Failed with MySQL on "<<mdbhost<<":"<<mdbPort<<endl;
      wm<<" Returned error =";
      wm<<mysql_error(&mData)<<".  Will re-try with timeout set at \n==> ";
      wm<<timeOutConnect<<" seconds <==";
      StDbManager::Instance()->printInfo((wm.str()).c_str(),dbMConnect,__LINE__,__CLASS__,__METHOD__); 
    }
  }      

  if(connected){
    if(mdbServerVersion) delete [] mdbServerVersion;
    mdbServerVersion=new char[strlen(mData.server_version)+1];
    strcpy(mdbServerVersion,mData.server_version);
  }

  return connected;

#undef __METHOD__
}

//////////////////////////////////////////////////////////////////////// 
bool MysqlDb::Connect(const char *aHost, const char *aUser, const char *aPasswd,  const char *aDb, const int aPort){
#define __METHOD__ "Connect(host,user,pw,database,port)"

  if(mdbhost) delete [] mdbhost;
  mdbhost  = new char[strlen(aHost)+1];   strcpy(mdbhost,aHost);
  if(aUser){
   if(mdbuser) delete [] mdbuser;
   mdbuser  = new char[strlen(aUser)+1];   strcpy(mdbuser,aUser);
   }
  if(aPasswd){
  if(mdbpw) delete [] mdbpw;
    mdbpw    = new char[strlen(aPasswd)+1]; strcpy(mdbpw,aPasswd);
  }
  mdbPort  = aPort;

  if(mdbName) {
    delete [] mdbName;
    mdbName=NULL;
  }
  if(aDb && (strcmp(aDb," ")!=0)){
    mdbName  = new char[strlen(aDb)+1];     
    strcpy(mdbName,aDb);
  }

  bool tRetVal = false;
  double t0=mqueryLog.wallTime();
  if(mlogTime)mconnectLog.start();
  if (!mysql_init(&mData))
    return (bool) StDbManager::Instance()->printInfo("Mysql Init Error=",mysql_error(&mData),dbMErr,__LINE__,__CLASS__,__METHOD__);

  // char *connString; 
  ostringstream cs;
  if(reConnect()){
    //  if(mysql_real_connect(&mData,aHost,aUser,aPasswd,bDb,aPort,NULL,0)){ 
       t0=mqueryLog.wallTime()-t0;
       cs<< "Server Connecting:"; if(mdbName)cs<<" DB=" << mdbName ;
       cs<< "  Host=" << aHost <<":"<<aPort<<endl;
       cs<< " --> Connection Time="<<t0<<" sec   ";
       if(mdbServerVersion)cs<<" MysqlVersion="<<mdbServerVersion;
      
        StDbManager::Instance()->printInfo((cs.str()).c_str(),dbMConnect,__LINE__,__CLASS__,__METHOD__);
      tRetVal=true;
  } else {
      cs << "Making Connection to DataBase = " << aDb;
      cs << " On Host = " << aHost <<":"<<aPort;
      cs << " MySQL returned error " << mysql_error(&mData);
      StDbManager::Instance()->printInfo((cs.str()).c_str(),dbMConnect,__LINE__,__CLASS__,__METHOD__);
  }

  if(mlogTime)mconnectLog.end();
  mhasConnected = tRetVal;
  return tRetVal;
#undef __METHOD__
}

////////////////////////////////////////////////////////////////////////

char* MysqlDb::printQuery(){ return mQueryLast; };

////////////////////////////////////////////////////////////////////////
void MysqlDb::RazQuery() {
  
  if (mQueryLast){
    delete [] mQueryLast;
    mQueryLast=0;
  }
  if(mhasBinaryQuery){
     mQueryLast = new char[strlen(binaryMessage)+1];
     strcpy(mQueryLast,binaryMessage);
     if(mQuery)delete [] mQuery;
  } else {
     if(mQuery){
       mQueryLast=new char[strlen(mQuery)+1];
       strcpy(mQueryLast,mQuery);
       delete [] mQuery;
     }
  }

  mQuery = 0;//new char[1];
  mQueryLen=0;
  //  strcpy(mQuery,"");

  mhasBinaryQuery=false;
}

bool MysqlDb::checkForTable(const char* tableName){

   mRes->Release();
   mRes->mRes=mysql_list_tables(&mData,tableName);
   if(mRes->mRes==NULL) return false;
   mRes->Release();

return true;
};

////////////////////////////////////////////////////////////////////////

bool MysqlDb::ExecQuery(){
#define __METHOD__ "ExecQuery()"

mqueryState=false;

//cout<<mQuery<<endl;

if(mlogTime)mqueryLog.start();

int status=mysql_real_query(&mData,mQuery,mQueryLen);

  if( (status!=0) && ( mysql_errno(&mData)==CR_SERVER_GONE_ERROR || mysql_errno(&mData)==CR_SERVER_LOST ) ){
       StDbManager::Instance()->printInfo(mysql_error(&mData)," Lost server, will try to reconnect",dbMDebug,__LINE__,__CLASS__,__METHOD__); 
       if(reConnect())status=mysql_real_query(&mData,mQuery,mQueryLen);       
  }

  if(status!=0)
    return StDbManager::Instance()->printInfo(" Query Failed ",mysql_error(&mData),dbMErr,__LINE__,__CLASS__,__METHOD__);

  if(mlogTime)mqueryLog.end();
  mRes->Release();

  if(mlogTime)msocketLog.start();
  mRes->mRes=mysql_store_result(&mData);
  if(mlogTime)msocketLog.end();

  return mqueryState=true;
}

////////////////////////////////////////////////////////////////////////

MysqlDb &MysqlDb::operator<<( const char *aQuery){ 
   
  if (strcmp(aQuery,";")==0){
    ExecQuery();
    RazQuery();
  } else {

    if(!mQuery){
      mQueryLen=strlen(aQuery);
      mQuery = new char[mQueryLen+1];
      strcpy(mQuery,aQuery);
    } else {
     char* tQuery = new char[strlen(mQuery)+1];
     strcpy(tQuery,mQuery);    
     delete [] mQuery;
     mQuery = new char[mQueryLen+strlen(aQuery)+1];
     memcpy(mQuery,tQuery,mQueryLen);
     strcpy(&mQuery[mQueryLen],aQuery);    
     delete [] tQuery;
     mQueryLen=mQueryLen+strlen(aQuery);
    }

  }

  return *this;
}

////////////////////////////////////////////////////////////////////////

MysqlDb &MysqlDb::operator<<( const MysqlBin *aBin ){

  mhasBinaryQuery=true;

  char *tQuery = new char[mQueryLen+aBin->mLen+1];
  memcpy(tQuery,mQuery,mQueryLen);
  memcpy(&tQuery[mQueryLen],aBin->mBinData,aBin->mLen+1); // mBinData included null terminator
  //  tQuery[mQueryLen+aBin->mLen]='\0';
  if(mQuery)delete [] mQuery;
  mQuery=tQuery;
  mQueryLen=mQueryLen+aBin->mLen;  // always not include null terminator
 
  return *this;
};

////////////////////////////////////////////////////////////////////////

bool MysqlDb::InputStart(const char* table,StDbBuffer *aBuff, const char* colList, int nRows,bool& hasBinary){

  bool tRetVal=false;
  if(!table || !aBuff || !colList) return tRetVal;

  bool change=aBuff->IsClientMode();
  if(change) aBuff->SetStorageMode();

  *this << "select * from " << table << " where null"<< endsql;
  *this << "insert delayed into " << table << " ("<<colList<<") VALUES(";
  int i;

    char* tmpString=new char[strlen(colList)+1];
    strcpy(tmpString,colList);
    char *ptr1,*ptr2;
    jfields=0;    
    bool done = false;
    ptr1=tmpString;

    while(!done){
      if((ptr2=strstr(ptr1,","))){
	*ptr2='\0';
      } else {
        done=true;
      }
      if(*ptr1==' ')ptr1++;
      if(cnames[jfields]) delete [] cnames[jfields];
      cnames[jfields]=new char[strlen(ptr1)+1];
      strcpy(cnames[jfields],ptr1);
      if(!done){
	ptr1=ptr2+1;
        *ptr2=',';
      }
      jfields++;
    }
    delete [] tmpString;
    int nfields=NbFields();
    int fcount=0;
    hasBinary=false;

    for(int k=0;k<jfields;k++){
      for(i=0;i<nfields;i++)
        if(strcmp(mRes->mRes->fields[i].name,cnames[k])==0)break;

      if(i==nfields)continue;
      fcount++;
      isBlob[k]=( (IS_BLOB(mRes->mRes->fields[i].flags)) || 
                  (mRes->mRes->fields[i].type ==254) );
      isBinary[k]= (mRes->mRes->fields[i].flags&BINARY_FLAG);
      isSpecialType[k]=(mRes->mRes->fields[i].type ==254);

      if(isBinary[k])hasBinary=true;
    }
    if(fcount!=jfields) done=false;

    return done;

};

bool MysqlDb::InputRow(StDbBuffer* aBuff, int row){

  char* tVal;tVal=0;
  int len; 
  aBuff->SetStorageMode();
      if(row>0)*this<<"),(";
      int k;
       for(k=0;k<jfields;k++){
         if(k!=0)*this<<",";
	 if(isBlob[k]){
	   if(isBinary[k]){
	     if(!aBuff->ReadArray(tVal,len,cnames[k]))break;
	     *this<<"'"<<Binary(len,(float*)tVal)<<"'";
	   } else if(isSpecialType[k]) {
             if(!aBuff->ReadScalar(tVal,cnames[k]))break;
             *this<<"'"<<tVal<<"'";
	   } else {
             char** tVal2=0;
	     if(!aBuff->ReadArray(tVal2,len,cnames[k]))break;
	     tVal=CodeStrArray(tVal2,len);
             for(int jj=0;jj<len;jj++)if(tVal2[jj])delete []tVal2[jj];
             *this<<"'"<<tVal<<"'";
	   }
	 } else {
	   if(!aBuff->ReadScalar(tVal,cnames[k])) break;
           *this<<"'"<<tVal<<"'";
 	 }
       }

       aBuff->SetClientMode();

    if(k!=jfields){
       RazQuery();
       return false;
    }

    return true;
}

bool MysqlDb::InputEnd(){

bool tRetVal=false;
*this<<")"<<endsql;
if(mqueryState)tRetVal=true;
  return tRetVal;

};
////////////////////////////////////////////////////////////////////////

bool MysqlDb::Input(const char *table,StDbBuffer *aBuff){
  bool tRetVal=false;
  bool change=aBuff->IsClientMode();
  if (change) aBuff->SetStorageMode();
  aBuff->SetStorageMode();
  if (aBuff) {
    *this << "select * from " << table << " where null"<< endsql;
    *this << "insert into " << table << " set ";
    bool tFirst=true;
    char* tVal;tVal=0;
    int len; 
    unsigned i;

 for (i=0;i<NbFields();i++) {
    if ((IS_BLOB(mRes->mRes->fields[i].flags) ) || 
        mRes->mRes->fields[i].type ==254) {
      if (mRes->mRes->fields[i].flags&BINARY_FLAG) { // Binary
	    if (aBuff->ReadArray(tVal,len,mRes->mRes->fields[i].name)){
	       if (tFirst) {
               tFirst=false;
           } else {
               *this << ",";
           };
	         *this << mRes->mRes->fields[i].name << "='" << Binary(len,(float*)tVal)<<"'";
	   };
	 }else{  // text types
       if(mRes->mRes->fields[i].type==254){
 	     if (aBuff->ReadScalar(tVal,mRes->mRes->fields[i].name)) {
	       if (tFirst) {
               tFirst=false;
           } else {
               *this << ",";
           };
	       *this << mRes->mRes->fields[i].name << "='" << tVal << "'";
	     };
       } else {
	     char** tVal2=0;
	     if (aBuff->ReadArray(tVal2,len,mRes->mRes->fields[i].name)){
	       tVal=CodeStrArray(tVal2,len);
	       int j;for (j=0;j<len;j++) {if (tVal2[j]) delete [] tVal2[j];};
	       delete [] tVal2;
	       if (tFirst) {
               tFirst=false;
           } else {
               *this << ",";
           };
	       *this << mRes->mRes->fields[i].name << "='" << tVal<<"'";
         }
       };
	};
  } else {  // not binary nor text
	      if (aBuff->ReadScalar(tVal,mRes->mRes->fields[i].name)) {
	         if (tFirst) {
                 tFirst=false;
             } else {
                 *this << ",";
             };
	         *this << mRes->mRes->fields[i].name << "='" << tVal << "'";
	      };
      };
      if (tVal) delete [] tVal;tVal=0;
    };
    if (tFirst) { 
      RazQuery();
    } else {
      *this << endsql;
      if(mqueryState)tRetVal=true;
    };
  };
  //  if (!tRetVal) cout << "insert Failed"<< endl;
  if (change) aBuff->SetClientMode();
  aBuff->SetClientMode();
  return tRetVal;
}; 

////////////////////////////////////////////////////////////////////////

bool  MysqlDb::Output(StDbBuffer *aBuff){

  if(mlogTime)msocketLog.start();
  MYSQL_ROW tRow=mysql_fetch_row(mRes->mRes);
  if(!tRow) return false;
  unsigned long * lengths=mysql_fetch_lengths(mRes->mRes);
  unsigned tNbFields=NbFields();
  if(mlogTime)msocketLog.end();
  int i;
  bool tRetVal=false;
  bool change=aBuff->IsClientMode();
  if (change) aBuff->SetStorageMode();
  aBuff->SetStorageMode();

  for (i=0;i<(int)tNbFields;i++){
    if(tRow[i]){

      //      cout <<mRes->mRes->fields[i].name<<" = "<< tRow[i] << endl;
      if (IS_BLOB(mRes->mRes->fields[i].flags)) {
	    if (mRes->mRes->fields[i].flags&BINARY_FLAG) {
	       aBuff->WriteArray((char*)tRow[i],lengths[i],mRes->mRes->fields[i].name);
	    }else {
	       char** tStrPtr;
	       int len;
	       tStrPtr=DecodeStrArray((char*)tRow[i],len);
	       aBuff->WriteArray(tStrPtr,len,mRes->mRes->fields[i].name);
           for(int k=0;k<len;k++)delete [] tStrPtr[k];
           delete [] tStrPtr;

           // something of a cludge: "," are used in "Decode.."  
           // as array delimeters. But we also want to have
           // char arrays for comment structures - so write'em both
           // to the buffer. 

           //char commentName[1024];
           ostringstream cn;
           cn<<mRes->mRes->fields[i].name<<".text";
	       aBuff->WriteScalar((char*)tRow[i],(cn.str()).c_str());
	    };
      } else {
	       aBuff->WriteScalar((char*)tRow[i],mRes->mRes->fields[i].name);
      };
    }

    /*
       else {

      ostringstream nd;
      nd<<"null data returned from table = ";
      nd<<mRes->mRes->fields[i].table<<" column="<<mRes->mRes->fields[i].name;
      StDbManager::Instance()->printInfo((nd.str()).c_str(),dbMWarn,__LINE__,"MysqlDb","Output(StDbBuffer* b)");
    }
    */

    tRetVal=true;
  };
  if (change) aBuff->SetClientMode();
  aBuff->SetClientMode();
  return tRetVal;
}

////////////////////////////////////////////////////////////////////////

char** MysqlDb::DecodeStrArray(char* strinput , int &aLen){
  
  if(!strinput){ // shouldn't happen - should have checked before here
    //   cout<< "null input string from mysql " << endl;
    char** tmparr = new char*[1];
    aLen = 1;
    *tmparr = new char[2];
    strcpy(*tmparr,"0");
    return tmparr;
  }

  char* tPnt=strinput; 
  aLen=0;
  while (tPnt&&aLen<100) { // 100 is a limit on # comma separated values
    tPnt=strpbrk( tPnt,"\\,");
    if (tPnt!=0){
      if (*tPnt==',') {
	aLen++;tPnt++;
      } else { //(*tPnt=='\\') 
	tPnt++;tPnt++;};
    };
  };
  aLen++;
  char** strarr=new char*[aLen];
  tPnt=strinput;
  char* tPntNew=tPnt;
  char *tBuff=new char[strlen(strinput)+1];
  char *tBuffInd=tBuff;
  int tCount=0;
  while (tPntNew) {
    tPntNew=strpbrk( tPnt,"\\,");
    if ((tPntNew==0)||*tPntNew==',') {
      if (tPntNew==0) {
	strcpy(tBuffInd,tPnt);
	
      } else {
	strncpy(tBuffInd,tPnt,tPntNew-tPnt);
	*(tBuffInd+(tPntNew-tPnt))='\0';
      }
      strarr[tCount]=new char[strlen(tBuff)+1];
      strcpy(strarr[tCount],tBuff);
      tBuffInd=tBuff;
      tPnt=tPntNew+1;
      tCount++;
    } else { //(*tPntNew=='\\') 
      strncpy(tBuffInd,tPnt,tPntNew-tPnt);
      tBuffInd=tBuffInd+(tPntNew-tPnt);
      tPntNew++;
      if(*tPntNew=='\\'||*tPntNew==',') {
	*tBuffInd=*tPntNew;
	tBuffInd++;
      };
      *(tBuffInd)='\0';
      tPnt=tPntNew+1;
    };
  };
  delete [] tBuff;
  return strarr;
}

////////////////////////////////////////////////////////////////////////

char* MysqlDb::CodeStrArray(char** strarr , int aLen){
  int tMaxLen=0;
  int i;
  for (i=0;i<aLen;i++) {
    if (strarr[i]) tMaxLen=tMaxLen+strlen(strarr[i])*2;
    tMaxLen++;
  };
  char* tTempVal=new char[tMaxLen+1];
  char* tRead;
  char* tWrite=tTempVal;
  for (i=0;i<aLen;i++) {
    if (strarr[i]){
      int j;
      tRead=strarr[i];
      for (j=0;j<(int)strlen(strarr[i]);j++) {
	if (*tRead=='\\'||*tRead==',') {
	  *tWrite='\\';
	  tWrite++;
	};
	*tWrite=*tRead;
	tWrite++;
	tRead++;
      };
    };
    *tWrite=',';
    tWrite++;
  };
  tWrite--;
  *tWrite='\0';
  char *tRetVal=new char[strlen(tTempVal)+1];
  strcpy(tRetVal,tTempVal);
  if (tTempVal) delete [] tTempVal;
  return tRetVal;
};
  
/////////////////////////////////////////////////////

bool
MysqlDb::setDefaultDb(const char* dbName){

  if(!dbName || strlen(dbName)==0)return false;
if(mdbName) delete [] mdbName;
mdbName=new char[strlen(dbName)+1];
strcpy(mdbName,dbName);
if(strcmp(dbName," ")==0)return true;

bool tOk=false;
unsigned int mysqlError;

 if(mysql_select_db(&mData,dbName)){

    mysqlError = mysql_errno(&mData);
    if(mysqlError==CR_SERVER_GONE_ERROR || mysqlError==CR_SERVER_LOST){
       reConnect();  
       if(mysql_select_db(&mData,dbName)){
         cerr<< "Error selecting database=" << dbName << endl;
         tOk=false;
       } else {
         tOk=true;
       }
    } else {
       cerr<< "Error selecting database=" << dbName << endl;
       tOk=false;
    }
 } else {
   tOk=true;
 }

return tOk;
}
  
////////////////////////////////////////////////////////////////////////

MysqlBin *Binary(const unsigned long int aLen,const float *aBin){
  
  MysqlBin *tBin=new MysqlBin;
  char *tString=new char[2*aLen+1];
  unsigned long int tNewLen=mysql_escape_string(tString,(char*) aBin,aLen);
  tBin->Input(tNewLen,tString);
  return tBin;
}

#undef __CLASS__


