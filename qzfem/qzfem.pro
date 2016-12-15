#-------------------------------------------------
#
# Project created by QtCreator 2014-07-30T09:11:00
#
#-------------------------------------------------

QT       -= gui

TARGET = qzfem
TEMPLATE = lib
CONFIG += staticlib

QMAKE_CXXFLAGS += -std=c++0x

# OpenMP section ####################################
QMAKE_CXXFLAGS += -DWITH_OPENMP # global definition for macro
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
#####################################################


SOURCES += hexahedralfem.cpp \
    mechanicalparameters3d.cpp \
    plasticfem.cpp \
    forcecondition3d.cpp \
    femforce.cpp \
    gaussquadrature.cpp \
    fem.cpp \
    planestressstrain.cpp \
    mindlinplatebending.cpp \
    fem2d.cpp \
    mindlinshellbending.cpp \
    mindlinplatelaminated.cpp

HEADERS += hexahedralfem.h \
    globalmatrix.h \
    femcondition.h \
    femcondition1d.h \
    femcondition2d.h \
    femcondition3d.h \
    mechanicalparameters3d.h \
    plasticfem.h \
    forcecondition3d.h \
    femforce.h \
    gaussquadrature.h \
    fem.h \
    planestressstrain.h \
    mindlinplatebending.h \
    fem2d.h \
    mindlinshellbending.h \
    mindlinplatelaminated.h
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

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../qzio/release/ -lqzio
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../qzio/debug/ -lqzio
else:unix: LIBS += -L$$OUT_PWD/../qzio/ -lqzio

INCLUDEPATH += $$PWD/../qzio
DEPENDPATH += $$PWD/../qzio

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzio/release/libqzio.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzio/debug/libqzio.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzio/release/qzio.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzio/debug/qzio.lib
else:unix: PRE_TARGETDEPS += $$OUT_PWD/../qzio/libqzio.a
