#!/bin/bash

set -ev

pre-commit run --all-files || ( git status --short ; git diff ; echo "** The files above do not have the correct indentation. **" ; echo "To fix, see instructions at https://github.com/wannier-developers/wannier90/wiki/ContributorsGuide#automatic-pre-commits-spaces-indentation-" ; exit 1 )
