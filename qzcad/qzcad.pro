#-------------------------------------------------
#
# Project created by QtCreator 2014-01-03T14:50:34
#
#-------------------------------------------------

QT       += core gui opengl xml script

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = qzcad
TEMPLATE = app

# OpenMP section ####################################
QMAKE_CXXFLAGS += -DWITH_OPENMP # global definition for macro
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp
#####################################################


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
    rotationbodymeshdialog.cpp \
    exportmeshdialog.cpp \
    qtxtsender.cpp \
    qstdredirector.cpp \
    boundaryconditionswidget.cpp \
    elasticfemdialog.cpp \
    qtscriptfemcondition3d.cpp \
    qtscriptforcecondition3d.cpp \
    elasticconstatntswidget.cpp \
    namedfloatingvector.cpp \
    codeeditor.cpp \
    qzscriptengine.cpp \
    qquadrilateralmesh2d.cpp \
    qpoint2d.cpp \
    qpoint3d.cpp \
    qtrianglemesh2d.cpp

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
    rotationbodymeshdialog.h \
    exportmeshdialog.h \
    qtxtsender.h \
    qstdredirector.h \
    boundaryconditionswidget.h \
    elasticfemdialog.h \
    qtscriptfemcondition3d.h \
    qtscriptforcecondition3d.h \
    elasticconstatntswidget.h \
    namedfloatingvector.h \
    codeeditor.h \
    qzscriptengine.h \
    qquadrilateralmesh2d.h \
    qpoint2d.h \
    qpoint3d.h \
    qtrianglemesh2d.h

FORMS    += mainwindow.ui \
    glcontrolwidget.ui \
    path2deditor.ui \
    structuredisomesh2ddialog.ui \
    polygonalmodeldialog.ui \
    pointeditordialog.ui \
    baryquadsdialog.ui \
    rectmeshdialog.ui \
    rotationbodymeshdialog.ui \
    exportmeshdialog.ui \
    boundaryconditionswidget.ui \
    elasticfemdialog.ui \
    elasticconstatntswidget.ui

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

OTHER_FILES += \
    license.txt

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../qzfem/release/ -lqzfem
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../qzfem/debug/ -lqzfem
else:unix: LIBS += -L$$OUT_PWD/../qzfem/ -lqzfem

INCLUDEPATH += $$PWD/../qzfem
DEPENDPATH += $$PWD/../qzfem

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzfem/release/libqzfem.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzfem/debug/libqzfem.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzfem/release/qzfem.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$OUT_PWD/../qzfem/debug/qzfem.lib
else:unix: PRE_TARGETDEPS += $$OUT_PWD/../qzfem/libqzfem.a

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
