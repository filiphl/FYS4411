TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt

LIBS += -larmadillo -lblas -llapack

SOURCES += main.cpp \
    system.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    particle.cpp \
    WaveFunctions/wavefunction.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/randomuniform.cpp \
    Math/random.cpp \
    sampler.cpp \
    WaveFunctions/simplegaussian.cpp \
    WaveFunctions/interactinsimplegaussian.cpp \
    Hamiltonians/interactingharmonicoscillator.cpp \
    optimizer.cpp \
    WaveFunctions/heliumwavefunction.cpp \
    Hamiltonians/heliumhamiltonian.cpp \
    WaveFunctions/twobodyquantumdot.cpp \
    Hamiltonians/twobodyquantumdothamiltonian.cpp \
    slater.cpp \
    WaveFunctions/manybodyquantumdotwavefunction.cpp \
    Hamiltonians/manybodyquantumdothamiltonian.cpp

HEADERS += \
    system.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    particle.h \
    WaveFunctions/wavefunction.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    Math/random.h \
    sampler.h \
    WaveFunctions/simplegaussian.h \
    WaveFunctions/interactinsimplegaussian.h \
    Hamiltonians/interactingharmonicoscillator.h \
    optimizer.h \
    WaveFunctions/heliumwavefunction.h \
    Hamiltonians/heliumhamiltonian.h \
    WaveFunctions/twobodyquantumdot.h \
    Hamiltonians/twobodyquantumdothamiltonian.h \
    slater.h \
    WaveFunctions/manybodyquantumdotwavefunction.h \
    Hamiltonians/manybodyquantumdothamiltonian.h

#INCLUDEPATH += /usr/include/openmpi-x86_64
#QMAKE_CXX = /usr/bin/mpicxx #/usr/lib64/openmpi/bin/mpicxx #mpicxx
#QMAKE_CXX_RELEASE = $$QMAKE_CXX
#QMAKE_CXX_DEBUG = $$QMAKE_CXX
#QMAKE_LINK = $$QMAKE_CXX
#QMAKE_CC = /usr/bin/mpicxx #/usr/lib64/openmpi/bin/mpicc #mpicc

#QMAKE_CFLAGS += $$system(mpicc --showme:compile)
#QMAKE_LFLAGS += $$system(mpicxx --showme:link)
#QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
#QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

