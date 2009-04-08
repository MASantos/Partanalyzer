#!/bin/bash

MANDIR=..
MANH=./man.h
MAN=$MANDIR/MAN
EXE=partanalyzer

SRC=`grep VERSION= *.cc *.h| grep -v "^[c]*tags"| head -1 | awk -F ":" '{print $1}'`
test "x$SRC" == "x" && echo "$0 : ERROR : No version tag found!" && exit

VERSION="`grep VERSION= $SRC | sed 's@.*="\(alpha\s[()-.0-9a-zA-Z]*\)";@\1@g'`"

test "x$VERSION" == "x" && echo "$0 : ERROR : No version found : Check source file $SRC" && exit

cat $MANH | sed 's@_VERSION_@'"$VERSION"'@g' > $MAN

$EXE --help >> $MAN


