#!/bin/bash

## This will download QE and compile it
version=5.4.0
destfile=espresso-${version}.tar.gz
url=http://qe-forge.org/gf/download/frsrelease/211/968/$destfile
destdir=espresso

export ESPRESSO_ROOT="$(pwd)/$destdir"

if [ -d $destdir ]
then
    if [ -d $destdir ]
    then
	rm -rf ${destdir}.old
    fi
    mv $destdir ${destdir}.old
fi
mkdir -p $destdir

wget $url --output-document=$destfile && tar xzf espresso-${version}.tar.gz -C $destdir --strip-components=1 && cd $destdir && ./configure && make -j pw pp

