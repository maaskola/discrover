#!/bin/bash
# A script to generate a Debian .deb package.
# Run it from the root directory of the Discrover package.
# It uses the current branch to generate a source archive.

export VERSION=`git describe | sed -e "s/-/./g"`
echo $VERSION

ruby scripts/make-source-release.rb packages/deb/discrover_$VERSION.tar.gz

cp packages/deb/discrover_$VERSION.tar.gz packages/deb/discrover_$VERSION.orig.tar.gz

cd packages/deb

tar xf discrover_$VERSION.orig.tar.gz 

# the debian directory contains some necessary files
cp -r debian/ discrover-$VERSION/
cd discrover-$VERSION/
dch --create -v $VERSION -M --package discrover
debuild -us -uc
