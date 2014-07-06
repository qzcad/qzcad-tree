#-------------------------------------------------
#
# Project created by QtCreator 2014-01-03T14:50:34
#
#-------------------------------------------------

QT       += core gui opengl xml script

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = qzcad
TEMPLATE = app

#QMAKE_CXXFLAGS += -std=c++0x

SOURCES += main.cpp\
        mainwindow.cpp \
    glmeshpicture.cpp \
    glcontrolwidget.cpp \
    path2deditor.cpp \
    structuredisomesh2ddialog.cpp \
    polygonalmodeldialog.cpp \
    pointeditordialog.cpp \
    baryquadsdialog.cpp \
    highlighter.cpp \
    colorvaluemap.cpp \
    rectmeshdialog.cpp \
    rotationbodymeshdialog.cpp

HEADERS  += mainwindow.h \
    globalconsts.h \
    glmeshpicture.h \
    glcontrolwidget.h \
    path2deditor.h \
    structuredisomesh2ddialog.h \
    polygonalmodeldialog.h \
    pointeditordialog.h \
    baryquadsdialog.h \
    highlighter.h \
    colorvaluemap.h \
    rectmeshdialog.h \
    rotationbodymeshdialog.h

FORMS    += mainwindow.ui \
    glcontrolwidget.ui \
    path2deditor.ui \
    structuredisomesh2ddialog.ui \
    polygonalmodeldialog.ui \
    pointeditordialog.ui \
    baryquadsdialog.ui \
    rectmeshdialog.ui \
    rotationbodymeshdialog.ui

RESOURCES += \
    qzresources.qrc

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../mesh/release/ -lmesh
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../mesh/debug/ -lmesh
else:unix: LIBS += -L$$OUT_PWD/../mesh/ -lmesh

INCLUDEPATH += $$PWD/../mesh
DEPENDPATH += $$PWD/../mesh

win32:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../mesh/release/mesh.lib
else:win32:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../mesh/debug/mesh.lib
else:unix: PRE_TARGETDEPS += $$OUT_PWD/../mesh/libmesh.a
