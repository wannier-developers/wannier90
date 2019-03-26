#!/usr/bin/env python
import subprocess
import sys

# Run the 'fprettify command with the same command line parameters
new_command = ['fprettify'] + sys.argv[1:]
# We put everything on one stream otherwise we should use threads
# to check and avoid locks
proc = subprocess.Popen(new_command, 
    stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
# Get stdout and stderr
outs, _ = proc.communicate()
# I write anyway to stderr
if outs is not None:
    if hasattr(sys.stderr, 'buffer'): # for Py3
        sys.stderr.buffer.write(outs)
    else:
        sys.stderr.write(outs)

# Return non-zero if there is any output. Replicate output
# (Note that output is printed as it comes)
retcode = proc.returncode
if retcode == 0 and outs:
        retcode = 200
    
sys.exit(retcode)
