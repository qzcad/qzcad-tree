TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

# OpenMP section ####################################
QMAKE_CXXFLAGS += -DWITH_OPENMP # global definition for macro

msvc {
  QMAKE_CXXFLAGS += -openmp
}

unix {
  QMAKE_CXXFLAGS += -fopenmp
  QMAKE_LFLAGS += -fopenmp
}

win32-g++ {
  QMAKE_CXXFLAGS += -fopenmp
  QMAKE_LFLAGS += -fopenmp
}
#####################################################


SOURCES += main.cpp


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
