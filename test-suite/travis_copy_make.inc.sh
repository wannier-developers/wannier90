#!/bin/bash
#Stop the full script if one line crashes
set -e

if [ "$W90BINARYPARALLEL" == "true" ]
then
  cp config/make.inc.gfort.traviscimpi make.inc
else
 cp config/make.inc.gfort.travisci make.inc 
fi
