#! /opt/star/bin/perl -w
#
# 
#
# 
#
# L.Didenko
#
# dbSumDAQProd.pl
#
# Scanning FilesCatalog table to get production summary and put it to Web page
# 
################################################################################################

use CGI;

require "/afs/rhic/star/packages/DEV00/mgr/dbCpProdSetup.pl";

use File::Find;
use Class::Struct;

my $debugOn=0;

#&cgiSetup();

my @SetD = (
             "P00hg/2000/06",
             "P00hg/2000/07",
             "P00hg/2000/08",
);
my $prod_html_back = "/star/u2e/starreco/P00hg/summary/P00hg.summary\.html_back";
my $prod_html = "/star/u2e/starreco/P00hg/summary/P00hg.summary\.html";

struct JFileAttr => {
    daqRun      => '$',
     prSer      => '$',  
    smFile      => '$',
    Nevts       => '$',
    jbStat      => '$', 
    EvtSk       => '$',
};

struct FilAttr => {
       flName   => '$',
       DtSt     => '$', 
       fpath    => '$',
       Nevts    => '$',
       numRun   => '$',
}; 
 
my @prodSer = ("P00hg","P00he") ;

&beginHtml();


##### Find sets in DataSet table

 my %dstDEvts = ();
 my %dstHpEvts = ();
 my %daqHpEvts = ();
 my %prodRun = ();
 my %periodRun = ();
 my %dtSet = ();
 my @DRun;
 my $nRun = 0;
 my $myRun;
 my @jobSum_set;
 my $jobSum_no = 0;
 my @dirP;
 my $dirR;
 my $topHpss = "/home/starreco/reco";
 my $topDisk = "/star/rcf/data09/reco"; 


  my %RunHash = (
                    "2000/06" => "JUNE-2000",
                    "2000/07" => "JULY-2000",
                    "2000/08" => "AUGUST-2000",
                    "2000/09" => "SEPTEMBER-2000"                    
 ); 

my %daqHash = ();


&StDbProdConnect();
##### select Geant files from FileCatalog
my  $nmfile;
my  @hpssInFiles;


 $sql="SELECT DISTINCT runID FROM $FileCatalogT WHERE jobID like '%P00hg%' AND fName LIKE '%dst.root' ";

 $cursor =$dbh->prepare($sql)
   || die "Cannot prepare statement: $DBI::errstr\n";
 $cursor->execute;

 while(@fields = $cursor->fetchrow) {
   my $cols=$cursor->{NUM_OF_FIELDS};

   for($i=0;$i<$cols;$i++) {
     my $fvalue=$fields[$i];
     my $fname=$cursor->{NAME}->[$i];
     print "$fname = $fvalue\n" if $debugOn;
      $myRun = $fvalue  if($fname eq 'runID');
#   print "runID : ", $myRun ,"\n";
   }
  
   $DRun[$nRun] = $myRun;
   $prodRun{$myRun} = $prodSer[0]; 
   $nRun++;
   }


 $sql="SELECT prodSeries, sumFileName, NoEvents, jobStatus, NoEventSkip FROM $JobStatusT WHERE prodSeries = '$prodSer[0]' and jobStatus <> 'n/a' ";


   $cursor =$dbh->prepare($sql)
    || die "Cannot prepare statement: $DBI::errstr\n";
   $cursor->execute;
 
    while(@fields = $cursor->fetchrow) {
     my $cols=$cursor->{NUM_OF_FIELDS};
        $fObjAdr = \(JFileAttr->new());
 

    for($i=0;$i<$cols;$i++) {
     my $fvalue=$fields[$i];
       my $fname=$cursor->{NAME}->[$i];
#        print "$fname = $fvalue\n" ;

        ($$fObjAdr)->prSer($fvalue)    if( $fname eq 'prodSeries');
        ($$fObjAdr)->jbStat($fvalue)   if( $fname eq 'jobStatus'); 
        ($$fObjAdr)->smFile($fvalue)   if( $fname eq 'sumFileName'); 
        ($$fObjAdr)->Nevts($fvalue)    if( $fname eq 'NoEvents');
        ($$fObjAdr)->EvtSk($fvalue)    if( $fname eq 'NoEventSkip')
   }

       $jobSum_set[$jobSum_no] = $fObjAdr;
       $jobSum_no++; 
}


##### select DST files on HPSS from FileCatalog
 my $hdfile;
 my $dhRun;
 my $dqRun;
 my $dhSet;
 my @OnlFiles;
 my $nOnlFile = 0;
 my @hpssDstFiles;
  $nhpssDstFiles = 0;

  for ($ll=0; $ll<scalar(@DRun); $ll++) {

  $sql="SELECT runID, dataset, fName, path, Nevents  FROM $FileCatalogT WHERE runID = '$DRun[$ll]' AND fName LIKE '%dst.root' AND JobID LIKE '%$prodSer[0]%' AND hpss = 'Y'";
  $cursor =$dbh->prepare($sql)
    || die "Cannot prepare statement: $DBI::errstr\n";
  $cursor->execute;

  while(@fields = $cursor->fetchrow) {
     my $cols=$cursor->{NUM_OF_FIELDS};
     $fObjAdr = \(FilAttr->new());

     for($i=0;$i<$cols;$i++) {
       my $fvalue=$fields[$i];
        my $fname=$cursor->{NAME}->[$i];
      print "$fname = $fvalue\n" if $debugOn;
 
      ($$fObjAdr)->DtSt($fvalue)     if( $fname eq 'dataset'); 
      ($$fObjAdr)->flName($fvalue)   if( $fname eq 'fName');
      ($$fObjAdr)->fpath($fvalue)    if( $fname eq 'path'); 
      ($$fObjAdr)->Nevts($fvalue)    if( $fname eq 'Nevents');
      ($$fObjAdr)->numRun($fvalue)   if( $fname eq 'runID');
    }
  
    $hpssDstFiles[$nhpssDstFiles] = $fObjAdr;  
    $nhpssDstFiles++;
   
    }
  }

 foreach my $dsfile (@hpssDstFiles) {

    $dhfile = ($$dsfile)->flName;
    $dhpath = ($$dsfile)->fpath;
    $dhRun = ($$dsfile)->numRun;
    @dirP = split ("/", $dhpath);
    $dirR = $dirP[5] . "/" . $dirP[6];

  $periodRun{$dhRun} = $RunHash{$dirR};
  $dstHpEvts{$dhRun}  += ($$dsfile)->Nevts; 

  }

##### select daq files from FileCatalog
  my $dqfile;
  my $dqRun;

  for ($ll=0; $ll<scalar(@DRun); $ll++) {

  $sql="SELECT runID, fName,path, Nevents  FROM $FileCatalogT WHERE runID = '$DRun[$ll]' AND fName LIKE '%daq' AND hpss ='Y'";
  $cursor =$dbh->prepare($sql)
     || die "Cannot prepare statement: $DBI::errstr\n";
  $cursor->execute;

  while(@fields = $cursor->fetchrow) {
    my $cols=$cursor->{NUM_OF_FIELDS};
    $fObjAdr = \(FilAttr->new());

    for($i=0;$i<$cols;$i++) {
      my $fvalue=$fields[$i];
      my $fname=$cursor->{NAME}->[$i];
      print "$fname = $fvalue\n" if $debugOn;

      ($$fObjAdr)->flName($fvalue)   if( $fname eq 'fName');
      ($$fObjAdr)->fpath($fvalue)   if( $fname eq 'path');     
      ($$fObjAdr)->Nevts($fvalue)    if( $fname eq 'Nevents');
      ($$fObjAdr)->numRun($fvalue)   if( $fname eq 'runID');
    }
    $OnlFiles[$nOnlFile] = $fObjAdr;  
    $nOnlFile++;
   
    }
  }


 foreach my $onfile (@OnlFiles) {

    $dqfile = ($$onfile)->flName;
    $dqpath = ($$onfile)->fpath;
    $dqRun = ($$onfile)->numRun;
    @dirP = split ("/", $dqpath);
    $dirR = $dirP[5] . "/" . $dirP[6];
   $periodRun{$dqRun} = $RunHash{$dirR}; 
   $daqHpEvts{$dqRun}  += ($$onfile)->Nevts; 
#  print "Period of run :", $dirR," % ", $periodRun{$dqRun}, "\n";

 }

#####  select DST files on DISK from FileCatalog

 my $ddfile;
 my $ddSet;
 my @diskDstFiles;

  $ndiskDstFiles = 0;

 for ($kk=0; $kk<scalar(@DRun); $kk++) {

  $sql="SELECT runID, fName, dataset, path, Nevents FROM $FileCatalogT WHERE runID = '$DRun[$kk]' AND fName LIKE '%dst.root' AND jobID LIKE '%$prodSer[0]%' AND site = 'disk_rcf'";
  $cursor =$dbh->prepare($sql)
    || die "Cannot prepare statement: $DBI::errstr\n";
  $cursor->execute;

  while(@fields = $cursor->fetchrow) {
     my $cols=$cursor->{NUM_OF_FIELDS};
     $fObjAdr = \(FilAttr->new());

    for($i=0;$i<$cols;$i++) {
       my $fvalue=$fields[$i];
       my $fname=$cursor->{NAME}->[$i];
       print "$fname = $fvalue\n" if $debugOn;

      ($$fObjAdr)->DtSt($fvalue)     if( $fname eq 'dataset');
      ($$fObjAdr)->flName($fvalue)    if( $fname eq 'fName');
      ($$fObjAdr)->fpath($fvalue)    if( $fname eq 'path');  
      ($$fObjAdr)->Nevts($fvalue)     if( $fname eq 'Nevents');
      ($$fObjAdr)->numRun($fvalue)    if( $fname eq 'runID');
    }
  
    $diskDstFiles[$ndiskDstFiles] = $fObjAdr;  
    $ndiskDstFiles++;
   }
  }
 
foreach my $ddfile (@diskDstFiles) {

   $dhSet  = ($$ddfile)->DtSt; 
   $dhfile = ($$ddfile)->flName;
   $dhpath = ($$ddfile)->fpath;
   $ddSet = ($$ddfile)->numRun;
    $dtSet{$ddSet} = $dhSet;
   $dstDEvts{$ddSet}  += ($$ddfile)->Nevts;   

  }


#initialize for total amount

my $TdstDEvt  = 0; 
my $TdstHEvt  = 0;
my $TdaqHEvt  = 0;


    foreach my $runD (@DRun) {
      
#        $daqHpEvts{$runD} -= 1;
        if (! defined $dstHpEvts{$runD}) {$dstHpEvts{$runD} = 0 };
        if (! defined $dstDEvts{$runD}) {$dstDEvts{$runD} = 0 };
        if (! defined $daqHpEvts{$runD}) {$daqHpEvts{$runD} = 0 };
        if (! defined $dtSet{$runD}) {$dtSet{$runD} = 'n/a'};
        $TdstHEvt    +=  $dstHpEvts{$runD}; 
        $TdstDEvt    +=  $dstDEvts{$runD}; 
        $TdaqHEvt    +=  $daqHpEvts{$runD}; 
#  print "Run ID", $runD, $periodRun{$runD}, "\n";

   print HTML "<TR ALIGN=CENTER VALIGN=CENTER>\n";
   print HTML "<TD><a href=\"http://duvall.star.bnl.gov/devcgi/dbP00hfDAQRetrv.pl?run=$runD\">$runD</TD>\n"; 
   print HTML "<td>$periodRun{$runD}</td><td>$prodRun{$runD}</td><td> $dtSet{$runD} </td><td>$daqHpEvts{$runD}</td><td>$dstHpEvts{$runD}</td><td>$dstDEvts{$runD}</td></tr>\n"; 

}

#  print total amount
   print HTML "<TR ALIGN=CENTER VALIGN=CENTER>\n"; 
   print HTML "<td>Total</td><td>June-August</td><td>All</td><td> All</td><td> $TdaqHEvt</td><td>$TdstHEvt</td><td>$TdstDEvt </td></tr>\n"; 

#####


##### finished with database
 &StDbProdDisconnect();

 &endHtml();

 `cp $prod_html_back $prod_html` ;

######################
 sub beginHtml {

   open (HTML,">$prod_html_back") or die "can't write to $prod_html_back ";
   print HTML "<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML//EN\">\n";
   print HTML "<html>\n";
   print HTML "  <head>\n";
   print HTML "          <title> Production summary </title>\n";
   print HTML "  </head>\n";
   print HTML "  <body BGCOLOR=\"#ccffff\"> \n";
   print HTML "      <h1>Production Summary for Summer 2000 </h1>\n";
   print HTML "<TABLE BORDER=5 CELLSPACING=1 CELLPADDING=2 >\n";
   print HTML "<TR>\n";
   print HTML "<TR ALIGN=CENTER VALIGN=CENTER>\n";
   print HTML "<TD WIDTH=\"10%\" HEIGHT=110><B>Run Number</B></TD>\n";
   print HTML "<TD WIDTH=\"20%\" HEIGHT=110><B>Month/Year</B></TD>\n";
   print HTML "<TD WIDTH=\"20%\" HEIGHT=110><B>Production period</B></TD>\n";
   print HTML "<TD WIDTH=\"20%\" HEIGHT=110><B>Dataset</B></TD>\n";
   print HTML "<TD WIDTH=\"10%\" HEIGHT=110><B>Number of Events<br> in DAQ files </B></TD>\n"; 
   print HTML "<TD WIDTH=\"10%\" HEIGHT=110><B>Number of Events<br>in DST on HPSS </B></TD>\n";
   print HTML "<TD WIDTH=\"10%\" HEIGHT=110><B>Number of Events<br>in DST on disk </B></TD>\n";
   print HTML "</TR>\n";
}

#####################
sub endHtml {
  
  my $Date = `date`;
  
  print HTML "</TABLE>\n";
  print HTML "      <h5>\n";
  print HTML "      <address><a href=\"mailto:didenko\@bnl.gov\">Lidia Didenko</a></address>\n";
  print HTML "<!-- Created: Mon Nov 29  05:29:25 MET 1999 -->\n";
  print HTML "<!-- hhmts start -->\n";
  print HTML "Last modified: $Date\n";
  print HTML "<!-- hhmts end -->\n";
  print HTML "  </body>\n";
  print HTML "</html>\n";
  close (HTML);
}












