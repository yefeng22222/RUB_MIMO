#-------------------------------------------------
#
# Project created by QtCreator 2015-09-16T14:51:25
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = OFDM
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    qcustomplot.cpp \
    figure.cpp \
    usrp_device.cpp \
    jsoncpp.cpp

HEADERS  += mainwindow.h \
    qcustomplot.h \
    types.h \
    figure.h \
    usrp_device.h \
    json/json.h \
    json/json-forwards.h

RESOURCES += \
    ressources.qrc

FORMS += \
    mainwindow.ui \
    no_dev.ui


unix: LIBS += -L/target/lib/ -luhd

INCLUDEPATH += /target/include
DEPENDPATH += /target/include
