#!/bin/csh -f
#       $Id$
#
#  D. Sandwell FEB 4 2010
#  M. Wei MAY 4 2010 - ENVISAT
#  E. Lindsey August 2017 - ALOS2 and TSX
# 
# Align a slave image to a master image and check results
#
#
alias rm 'rm -f'
unset noclobber
#
# check the number of arguments 
# 
  if ($#argv < 3) then 
    echo ""
    echo "Usage: align.csh SAT master_name slave_name [supermaster_name]"
    echo ""
    echo " The supermaster_namestem is required if this is secondary alignment."
    echo " SAT = ERS, ENVI, ALOS, ALOS2, TSX or generic SAT"
    echo ""
    echo "Example: align.csh ALOS IMG-HH-ALPSRP055750660-H1.0__A IMG-HH-ALPSRP049040660-H1.0__A "
    echo ""
    exit 1
  endif
  set SAT = $1
  if( ($SAT != ALOS2) && ($SAT != ALOS) && ($SAT != TSX) && ($SAT != ENVI) && ($SAT != ERS) && ($SAT != SAT)) then
    echo ""
    echo " SAT must be ERS, ENVI, ALOS, ALOS2, TSX, or generic SAT"
    echo ""
    exit 1
  endif
#
# ALOS2 is a special case
#
if ($SAT == ALOS2) then

# warning: supermaster is not being handled for ALOS2 at all
echo "align.csh: warning: ALOS-2 alignment ignores the supermaster"
echo ""

  set mode = $5
  if ($mode == scan) then

    cp $2.PRM $2.PRM0
    cp $3.PRM $3.PRM0

    # update slave.PRM 
    echo "update slave.PRM"
    set RSHIFT = `SAT_baseline $2.PRM $3.PRM | grep rshift | awk '{print $3}'`
    set ASHIFT = `SAT_baseline $2.PRM $3.PRM | grep ashift | awk '{print $3}'`
    update_PRM $3.PRM rshift $RSHIFT
    update_PRM $3.PRM ashift $ASHIFT

    echo "correlate master and slave to find offset parameters"
    xcorr $2.PRM $3.PRM -xsearch 64 -ysearch 64 -nx 32 -ny  128

    mv $3.PRM junk.PRM
    cp $2.PRM0 $2.PRM
    grep -v shift < junk.PRM > $3.PRM

    awk '{print $4}' < freq_xcorr.dat > tmp.dat
    set amedian = `sort -n tmp.dat | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'`
    set amax = `echo $amedian | awk '{print $1+3}'`
    set amin = `echo $amedian | awk '{print $1-3}'`
    awk '{if($4 > '$amin' && $4 < '$amax') print $0}' < freq_xcorr.dat > tmp2.dat

    fitoffset.csh 2 2 tmp2.dat 25 >> $3.PRM
  
  else

    xcorr $2.PRM $3.PRM -xsearch 64 -ysearch 64 -nx 32 -ny 64
    fitoffset.csh 2 2 freq_xcorr.dat 18 >> $3.PRM

  endif 

  resamp $2.PRM $3.PRM $3.PRMresamp $3.SLCresamp 4

  rm $3.SLC
  mv $3.SLCresamp $3.SLC
  cp $3.PRMresamp $3.PRM

  # end ALOS2 special case

else if ($SAT == TSX) then
  # TSX special case
  cp $3.PRM $3.PRM0
  SAT_baseline $2.PRM $3.PRM0 >> $3.PRM
  xcorr $2.PRM $3.PRM -xsearch 128 -ysearch 128
  fitoffset.csh 2 2 freq_xcorr.dat >> $3.PRM
  resamp $2.PRM $3.PRM $3.PRMresamp $3.SLCresamp 4
  rm $3.SLC
  mv $3.SLCresamp $3.SLC
  cp $3.PRMresamp $3.PRM

else
  
#
# focus the master if necessary
# Do it no matter what for now. Put SLC_file to PRM. Might not be necessary
#
  if(! -f $2.SLC) then
    echo "focussing master"
    sarp.csh $2.PRM 
  else
    update_PRM $2.PRM SLC_file $2.SLC
  endif
#
# focus the slave image
#
# check the range sampling rate 
# 
  set rng_samp_rate_m = `grep rng_samp_rate $2.PRM | awk 'NR == 1 {printf("%d", $3)}'`
  set rng_samp_rate_s = `grep rng_samp_rate $3.PRM | awk 'NR == 1 {printf("%d", $3)}'`
  if ($rng_samp_rate_m != $rng_samp_rate_s) then 
    echo "The range sampling rate for master and slave differ"
    echo "Need to run the interferogram in steps until process2pass.csh is fixed"
    exit 1
  endif 
  echo "align.csh"
  echo "focusing slave"
  sarp.csh $3.PRM 
#
# get the starting alignment parameters and run xcorr
#
  cp $2.PRM $2.PRM0
  cp $3.PRM $3.PRM0
  if($#argv == 4) then
    #set RSHIFT = `$1_baseline $4.PRM $3.PRM | grep rshift | awk '{print $3}'`
    #set ASHIFT = `$1_baseline $4.PRM $3.PRM | grep ashift | awk '{print $3}'`
    set RSHIFT = `SAT_baseline $4.PRM $3.PRM | grep rshift | awk '{print $3}'`
    set ASHIFT = `SAT_baseline $4.PRM $3.PRM | grep ashift | awk '{print $3}'`
#
#   use the PRF of the supermaster in the surrogate master
#
    set PRF = `grep PRF $4.PRM | awk '{print $3}'`
    update_PRM $2.PRM PRF $PRF
  else
    #set RSHIFT = `$1_baseline $2.PRM $3.PRM | grep rshift | awk '{print $3}'`
    #set ASHIFT = `$1_baseline $2.PRM $3.PRM | grep ashift | awk '{print $3}'`
    set RSHIFT = `SAT_baseline $2.PRM $3.PRM | grep rshift | awk '{print $3}'`
    set ASHIFT = `SAT_baseline $2.PRM $3.PRM | grep ashift | awk '{print $3}'`
  endif
  update_PRM $3.PRM rshift $RSHIFT
  update_PRM $3.PRM ashift $ASHIFT
  echo "align.csh"
  echo "correlate master and slave to find offset parameters"
  if( $SAT == "ERS") then
    xcorr $2.PRM $3.PRM -xsearch 128 -ysearch 128 -nx 20 -ny 50
  else
    xcorr $2.PRM $3.PRM -xsearch 128 -ysearch 256 -nx 20 -ny 50
  endif
#
  mv $3.SLC $3.SLC0
  mv $3.PRM junk.PRM
  cp $2.PRM0 $2.PRM
  grep -v shift < junk.PRM > $3.PRM
#
# put in the alignment parameters 
#
  fitoffset.csh 3 3 freq_xcorr.dat 18 >> $3.PRM
  mv freq_xcorr.dat xcorr_$2_$3.dat0
#
# refocus the second image
#
  echo "align.csh"
  echo "refocus slave"
  sarp.csh $3.PRM 
#
rm *SLC0
rm junk*

endif #end non-ALOS2 processing case
#done
