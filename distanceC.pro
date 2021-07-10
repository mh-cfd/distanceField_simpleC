TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE += -O3

QMAKE_LFLAGS += -O3

SOURCES += main.c

HEADERS += \
    calc_distance.h
