INCLUDEPATH += /c/PFM_ABEv7.0.0_Win64/include
LIBS += -L /c/PFM_ABEv7.0.0_Win64/lib -lnvutility -lgdal -lxml2 -lpoppler -lz -lm -liconv
DEFINES += NVWIN3X
CONFIG += console
CONFIG -= qt
QMAKE_LFLAGS += 
######################################################################
# Automatically generated by qmake (2.01a) Wed Jan 22 13:43:28 2020
######################################################################

TEMPLATE = app
TARGET = build_swbd
DEPENDPATH += .
INCLUDEPATH += .

# Input
HEADERS += version.h
SOURCES += main.c
