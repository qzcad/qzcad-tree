#-------------------------------------------------
#
# Project created by QtCreator 2015-01-02T12:02:31
#
#-------------------------------------------------

QT       -= core gui

TARGET = qzio
TEMPLATE = lib
CONFIG += staticlib

SOURCES += consoleprogress.cpp

HEADERS += consoleprogress.h
unix {
    target.path = /usr/lib
    INSTALLS += target
}
