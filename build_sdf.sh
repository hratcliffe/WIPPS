#!/bin/bash

export WKDIR=`pwd`
cd $1
make
cd $WKDIR
