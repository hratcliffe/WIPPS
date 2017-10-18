#!/bin/bash
#Get the current code version (if a git repo), or a hard-coded string version otherwise

#Hardcoded string if downloaded not cloned
HARDCODE='v1.1-00-000000'

GIT_VERSION=`git describe --always --tags | cut -c1-14`
# This encodes the git commit version into the source so we can write version number into our data files etc

if [ $? -ne 0 ]; then
  GIT_VERSION=$HARDCODE
fi

echo 'GIT_VERSION='$GIT_VERSION > `dirname $0`"/VERSION"
