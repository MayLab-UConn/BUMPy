#!/bin/bash

# unit tests are not ideally set up in terms of relative paths, so we'll cd into 
# directories and run tests from there

topdir=`pwd`
cd $topdir/unittests
echo "Running unit tests"
python3 -m unittest
if [ $? -eq 0 ]; then
    echo "Unit tests passed."
    echo "Running command line tests."
    cd $topdir/python_cli_tests
    python3 -m unittest
    if [ $? -eq 0 ]; then
        echo "Command line tests passed"
    else
	echo "Command line tests failed"
	cd $topdir
	exit 1
    fi
else
    echo "Unit tests failed"
    cd $topdir
    exit 1
fi
echo "All tests passed!"
cd $topdir
exit 0
