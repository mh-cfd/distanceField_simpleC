TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CFLAGS_RELEASE += -O3 -msse2 -march=native

QMAKE_CFLAGS += -O3 -msse2 -march=native

QMAKE_LFLAGS += -O3 -msse2 -march=native

SOURCES += main.c

HEADERS += \
    calc_distance.h
