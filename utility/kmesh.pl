#!/usr/bin/perl -w
#
# This file is distributed as part of the Wannier90 code and
# under the terms of the
#
#      GNU Lesser General Public License, version 2.1
#
# see file `LICENSE' in the root directory of the Wannier90
# distribution or https://www.gnu.org/licenses/lgpl-2.1.html
#
# The webpage of the Wannier90 code is www.wannier.org
#
# The Wannier90 code is hosted on GitHub
# https://github.com/wannier-developers/wannier90
#

$numargs = $#ARGV+1;

if (($numargs<3)||($numargs>4)) {
    print  "usage: n1 n2 n3 [wan]\n";
    print  "       n1  - divisions along 1st recip vector\n";
    print  "       n2  - divisions along 2nd recip vector\n";
    print  "       n3  - divisions along 3rd recip vector\n";
    print  "       wan - omit the kpoint weight (optional)\n";
    exit;
}

if ($ARGV[0]<=0) {
    print "n1 must be >0\n";
    exit;
}
if ($ARGV[1]<=0) {
    print "n2 must be >0\n";
    exit;
}
if ($ARGV[2]<=0) {
    print "n3 must be >0\n";
    exit;
}

$totpts=$ARGV[0]*$ARGV[1]*$ARGV[2];

if ($numargs==3) {
    print "K_POINTS crystal\n";
    print $totpts,"\n";
    for ($x=0; $x<$ARGV[0]; $x++) {
	for ($y=0; $y<$ARGV[1]; $y++) {
	    for ($z=0; $z<$ARGV[2]; $z++) {
		printf ("%12.8f%12.8f%12.8f%14.6e\n", $x/$ARGV[0],$y/$ARGV[1],$z/$ARGV[2],1/$totpts);
	    }
	}
    }
}


if ($numargs==4) {
    for ($x=0; $x<$ARGV[0]; $x++) {
	for ($y=0; $y<$ARGV[1]; $y++) {
	    for ($z=0; $z<$ARGV[2]; $z++) {
		printf ("%12.8f%12.8f%12.8f\n", $x/$ARGV[0],$y/$ARGV[1],$z/$ARGV[2]);
	    }
	}
    }
}


exit;
