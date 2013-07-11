#-------------------------------------------------
#
# Project created by QtCreator 2013-06-28T22:01:39
#
#-------------------------------------------------

QT       += core gui
QT       += opengl

TARGET = Stable_Fluids
TEMPLATE = app

CONFIG += c++11

INCLUDEPATH += C:\code\eigen-eigen-2249f9c22fe8


SOURCES += main.cpp\
        glwidget.cpp \
    simulation.cpp

HEADERS  += glwidget.h \
    vec.h \
    simulation.h \
    pcgsolver/sparse_matrix.h \
    pcgsolver/pcg_solver.h \
    pcgsolver/blas_wrapper.h \
    util.h \
    mac/types.h \
    mac/traits.h \
    mac/macgridfactory.h \
    mac/mac.h \
    mac/lerp.h \
    mac/grid.h
