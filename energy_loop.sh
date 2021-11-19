#!/bin/bash
set -o errexit

filename="mymac_WScFi.mac"
num_threads=1

particle="e-"

num_events=5000
energies=(1 2 5 10 20 30 40 50 60 70 80 90 100)
k=$1

mkdir -p ${particle}_QGSP
mkdir -p proc$1
cd proc$1

cp ../$filename .
sed -i "s/\/analysis\/setFileName .*/\/analysis\/setFileName out/" $filename
sed -i "s/\/gps\/particle .*/\/gps\/particle ${particle}/" $filename
sed -i "s/\/gps\/ene\/mono .*/\/gps\/ene\/mono ${energies[k]} GeV/" $filename
sed -i "s/\/run\/beamOn .*/\/run\/beamOn ${num_events}/" $filename
../build/ATHENA_Geometry -m mymac_WScFi.mac -t ${num_threads}
mv out.root ../${particle}_QGSP/${particle}_${energies[k]}GeV_homo.root

cd ..
rm -rf proc$1
