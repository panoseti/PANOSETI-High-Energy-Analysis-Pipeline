#!/usr/bin/env bash
RUNDIR=
IODIR=
OUTDIR=

THISDIR=$(pwd)

# remove old files
rm -f batch*.root
rm -f batch*.log
rm -f joblog.txt
rm -f $RUNDIR/DATbatch*
rm -f $IODIR/batch*.root

# run corsika and corsikaIOreader
parallel --verbose --jobs 16 --joblog joblog.txt ./corsika.sh {} {/.} ::: $THISDIR/batch*.inp

# copy files over
cd $OUTDIR
cp $IODIR/batch*.root .