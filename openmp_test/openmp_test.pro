TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

QMAKE_CXXFLAGS += -std=c++0x

# OpenMP section ####################################
QMAKE_CXXFLAGS += -DWITH_OPENMP # global definition for macro
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
#####################################################

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
