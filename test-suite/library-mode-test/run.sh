#!/bin/bash
make -j -C ../.. lib && make && ./test_library_serial.x
