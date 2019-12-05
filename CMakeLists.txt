cmake_minimum_required(VERSION 3.14)
project(Lab1Vorobiov C)

set(CMAKE_C_STANDARD 11)

add_executable(Lab1Vorobiov main.c
        bisection.c bisection.h
        discretenewton.c discretenewton.h
        newton.c newton.h
        differentialinterpolation.c differentialinterpolation.h
        cubicspline.c cubicspline.h
        bezier.c bezier.h
        matrixutils.c matrixutils.h
        householder.c householder.h
        quadraticregression.c quadraticregression.h
        bilagrange.c bilagrange.h
        simpleintegrals.c simpleintegrals.h)
