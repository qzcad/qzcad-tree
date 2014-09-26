#-------------------------------------------------
#
# Project created by QtCreator 2014-07-30T09:11:00
#
#-------------------------------------------------

QT       -= gui

TARGET = qzfem
TEMPLATE = lib
CONFIG += staticlib

SOURCES += hexahedralfem.cpp \
    mechanicalparameters3d.cpp

HEADERS += hexahedralfem.h \
    globalmatrix.h \
    femcondition.h \
    femcondition1d.h \
    femcondition2d.h \
    femcondition3d.h \
    mechanicalparameters3d.h
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

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../qzmatrix/release/ -lqzmatrix
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../qzmatrix/debug/ -lqzmatrix
else:unix: LIBS += -L$$OUT_PWD/../qzmatrix/ -lqzmatrix

INCLUDEPATH += $$PWD/../qzmatrix
DEPENDPATH += $$PWD/../qzmatrix

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzmatrix/release/libqzmatrix.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzmatrix/debug/libqzmatrix.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzmatrix/release/qzmatrix.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzmatrix/debug/qzmatrix.lib
else:unix: PRE_TARGETDEPS += $$OUT_PWD/../qzmatrix/libqzmatrix.a
