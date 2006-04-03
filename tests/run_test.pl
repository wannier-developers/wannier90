#!/usr/bin/perl
#
use Cwd;
my $cwd = cwd();
$wanex="../wannier90.x" ;
$wanex=$cwd."/".$wanex ;
print " Running test set of inputs\n";
@tests_found = <./test*>;
open (L,">wantest.log")
    || die " Can't open ".$cwd."\/wantest.log\n";

foreach $want (@tests_found){
    $tname = $want;
    $tname =~ s/\.\///; 
    chdir ("$want") or die "Can't cd to: $!\n" ; 
    print " Starting ".$tname."\n" ;
    print L "*" x 37 ." " . $tname ." ". "*" x 37 ."\n\n";
    open (D,"<des.dat")
        || print " Can't open ".$cwd."\/".$tname."\/des.dat\n";
    while (<D>) {
	print $_;
	print L " " x 15 . $_ ;
    }
    print L "\n";
    close (D);
    system("$wanex > /dev/null");
    $count = 0;
    unless (open(O,"<wannier.wout")) {
	$count=1;;
    }
    if ($count == 0) {
	chomp($last_line = $_) while(<O>);
	close (O);
	if ($last_line !~ /^*All done: wannier90 exiting*/) 
	{
	    $count =2;
	}
    }
    if ($count == 0) {
	print " Test ran to completion \n\n";
	print L "=" x36 ." Standard "."=" x36 . "\n\n";
	get_spread("stnd.wout");
	print L "=" x36 ." Current  "."=" x36 . "\n\n";
	get_spread("wannier.wout");
    } else {
	print " Test Failed. \n\n";
    }
    
    
    chdir ("../") ;
}
print " Tests Complete\n";
print " Examine the file ".$cwd."\/wantest.log\n";
close(L);

sub get_spread {
    my $file = shift;
    open (O,"<$file")
        || print " Can't open ".$cwd."\/".$tname."\/" . $file ."\n";
    while (<O>) {
	print L if /^*Final State*/ .. /^*Final Spread*/;
    }
    print L "\n\n";
    close(O);
}
