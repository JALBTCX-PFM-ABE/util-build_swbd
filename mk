#!/bin/bash

if [ ! $PFM_ABE_DEV ]; then

    export PFM_ABE_DEV=${1:-"/usr/local"}

fi

export PFM_BIN=$PFM_ABE_DEV/bin
export PFM_LIB=$PFM_ABE_DEV/lib
export PFM_INCLUDE=$PFM_ABE_DEV/include


CHECK_QT=`echo $QTDIR | grep "qt-3"`
if [ $CHECK_QT ] || [ !$QTDIR ]; then
    QTDIST=`ls ../../FOSS_libraries/qt-*.tar.gz | cut -d- -f5 | cut -dt -f1 | cut -d. --complement -f4`
    QT_TOP=Trolltech/Qt-$QTDIST
    QTDIR=$PFM_ABE_DEV/$QT_TOP
fi


SYS=`uname -s`


if [ $SYS = "Linux" ]; then
    DEFS="NVLinux"
    LIBRARIES="-L $PFM_LIB -lnvutility -lgdal -lxml2 -lpoppler -lz -lGLU -lm"
    export LD_LIBRARY_PATH=$PFM_LIB:$QTDIR/lib:$LD_LIBRARY_PATH
else
    DEFS="NVWIN3X"
    LIBRARIES="-L $PFM_LIB -lnvutility -lgdal -lxml2 -lpoppler -lz -lm -liconv"
    export QMAKESPEC=win32-g++
fi


# As of gcc 6 --enable-default-pie has been built in to the gcc compiler.
# We need to turn it off.

GVERSION=`gcc -dumpversion | cut -f 1 -d.`
MFLAGS=""
if [ $GVERSION -gt 5 ]; then
    MFLAGS=-no-pie
fi


#  Get the name from the directory name

NAME=`basename $PWD`


# Building the Makefile using qmake and adding extra includes, defines, and libs


rm -f $NAME.pro Makefile

$QTDIR/bin/qmake -project -o $NAME.tmp
cat >$NAME.pro <<EOF
INCLUDEPATH += $PFM_INCLUDE
LIBS += $LIBRARIES
DEFINES += $DEFS
CONFIG += console
CONFIG -= qt
QMAKE_LFLAGS += $MFLAGS
EOF

cat $NAME.tmp >>$NAME.pro
rm $NAME.tmp


$QTDIR/bin/qmake -o Makefile



if [ $SYS = "Linux" ]; then
    make
    if [ $? != 0 ];then
        exit -1
    fi
    chmod 755 $NAME
    mv $NAME $PFM_BIN
else
    if [ ! $WINMAKE ]; then
        WINMAKE=release
    fi
    make $WINMAKE
    if [ $? != 0 ];then
        exit -1
    fi
    chmod 755 $WINMAKE/$NAME.exe
    cp $WINMAKE/$NAME.exe $PFM_BIN
    rm $WINMAKE/$NAME.exe
fi


# Get rid of the Makefile so there is no confusion.  It will be generated again the next time we build.

rm Makefile
