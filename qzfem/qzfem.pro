#-------------------------------------------------
#
# Project created by QtCreator 2014-07-30T09:11:00
#
#-------------------------------------------------

QT       -= gui

TARGET = qzfem
TEMPLATE = lib
CONFIG += staticlib

SOURCES += hexahedralfem.cpp

HEADERS += hexahedralfem.h \
    globalmatrix.h \
    floatingvector.h \
    floatingmatrix.h
unix {
    target.path = /usr/lib
    INSTALLS += target
}

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../mesh/release/ -lmesh
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../mesh/debug/ -lmesh
else:unix: LIBS += -L$$OUT_PWD/../mesh/ -lmesh

INCLUDEPATH += $$PWD/../mesh
DEPENDPATH += $$PWD/../mesh

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../mesh/release/libmesh.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../mesh/debug/libmesh.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../mesh/release/mesh.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../mesh/debug/mesh.lib
else:unix: PRE_TARGETDEPS += $$OUT_PWD/../mesh/libmesh.a
