QT += core gui opengl
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

INCLUDEPATH += code
INCLUDEPATH += code/core
INCLUDEPATH += code/appRender
INCLUDEPATH += $$(GLEW_DIR)/include/GL

VPATH += code

SOURCES += \
    glacierterrain.cpp \
    appRender/main.cpp \
    appRender/mainwindow.cpp \
    appRender/glacierswidget.cpp \
    appRender/featurespainter.cpp \
    core/core.cpp \
    core/scalarfield2.cpp \
    core/vectorfield2.cpp \
    core/shader-utils.cpp

HEADERS += \
    glacierterrain.h \
    appRender/mainwindow.h \
    appRender/glacierswidget.h \
    appRender/featurespainter.h \
    core/core.h \
    core/scalarfield2.h \
    core/vectorfield2.h \
    core/shader-utils.h

FORMS += \
    appRender/mainwindow.ui

LIBS += -L$$(GLEW_DIR)/lib/Release/x64 -lglew32
LIBS += -lopengl32 -lglu32

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
