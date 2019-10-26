cmake_minimum_required(VERSION 3.14)
project(Lab1Vorobiov C)

set(CMAKE_C_STANDARD 11)

add_executable(Lab2 main.c
        bisection.c bisection.h
        discretenewton.c discretenewton.h
        newton.c newton.h
        differentialinterpolation.c differentialinterpolation.h)
