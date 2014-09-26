TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

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