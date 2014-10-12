#-------------------------------------------------
#
# Project created by QtCreator 2014-09-23T09:10:19
#
#-------------------------------------------------

QT       -= core gui

TARGET = qzmatrix
TEMPLATE = lib
CONFIG += staticlib

# OpenMP section ####################################
QMAKE_CXXFLAGS += -DWITH_OPENMP # global definition for macro
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
#####################################################

SOURCES += doublematrix.cpp \
    doublevector.cpp \
    mappeddoublematrix.cpp \
    rowdoublematrix.cpp

HEADERS += doublematrix.h \
    doublevector.h \
    mtx_types.h \
    mappeddoublematrix.h \
    rowdoublematrix.h
unix {
    target.path = /usr/lib
    INSTALLS += target
}
