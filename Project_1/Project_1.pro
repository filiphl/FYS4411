TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    solver.cpp \
    random.cpp

HEADERS += \
    solver.h \
    random.h

DISTFILES += \
    readdate.py

