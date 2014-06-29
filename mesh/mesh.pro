#-------------------------------------------------
#
# Project created by QtCreator 2014-01-03T14:51:08
#
#-------------------------------------------------

QT       -= core gui

TARGET = mesh
TEMPLATE = lib
CONFIG += staticlib
QMAKE_CXXFLAGS += -std=c++0x

SOURCES += mesh.cpp \
    point.cpp \
    element.cpp \
    point1d.cpp \
    point2d.cpp \
    quadrilateral.cpp \
    path2d.cpp \
    segment.cpp \
    mesh2d.cpp \
    quadrilateralmesh2d.cpp \
    quadrilateralunion2d.cpp \
    point3d.cpp \
    mesh3d.cpp \
    hexahedral.cpp \
    hexahedralmesh3d.cpp

HEADERS += mesh.h \
    point.h \
    floating.h \
    element.h \
    point1d.h \
    point2d.h \
    quadrilateral.h \
    pointpointer.h \
    elementpointer.h \
    meshpointer.h \
    path2d.h \
    segment.h \
    integer.h \
    mesh2d.h \
    quadrilateralmesh2d.h \
    adjacentset.h \
    node2d.h \
    nodetype.h \
    quadrilateralunion2d.h \
    point3d.h \
    integervector.h \
    node3d.h \
    mesh3d.h \
    hexahedral.h \
    hexahedralmesh3d.h
unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}
