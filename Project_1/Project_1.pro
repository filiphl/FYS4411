TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

SOURCES += main.cpp \
    solver.cpp \
    random.cpp

HEADERS += \
    solver.h \
    random.h

DISTFILES += \
    readdate.py

