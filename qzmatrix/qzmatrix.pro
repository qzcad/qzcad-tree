#-------------------------------------------------
#
# Project created by QtCreator 2014-09-23T09:10:19
#
#-------------------------------------------------

QT       -= core gui

TARGET = qzmatrix
TEMPLATE = lib
CONFIG += staticlib

SOURCES += doublematrix.cpp \
    doublevector.cpp

HEADERS += doublematrix.h \
    doublevector.h \
    mtx_types.h
unix {
    target.path = /usr/lib
    INSTALLS += target
}
