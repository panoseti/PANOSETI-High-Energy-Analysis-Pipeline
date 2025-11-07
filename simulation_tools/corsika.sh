#!/usr/bin/env bash
RUNDIR=
IODIR=
OUTDIR=

#corsika
cd $RUNDIR
./corsika77410Linux_QGSII_urqmd <$1 > "$OUTDIR/$2.log" &

wait $!

#corsikaIOreader
cd $IODIR
corsikaIOreader -cors "$RUNDIR/DAT$2.telescope" -histo "$2.root" -abs CORSIKA &

wait $!