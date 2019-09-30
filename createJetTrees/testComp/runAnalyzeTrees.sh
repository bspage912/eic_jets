#!/bin/sh

set -o noglob

./analyzeTest.exe 0 /gpfs02/eic/bpage/testWTA/jetTree.ep.20x250.1Mevents.RadCor=0.Q2=*.root test.hist.root
