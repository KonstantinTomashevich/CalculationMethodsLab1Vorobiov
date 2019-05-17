cmake_minimum_required(VERSION 3.14)
project(Lab2 C)

set(CMAKE_C_STANDARD 11)

add_executable(Lab2 main.c def.h
        matrixutils.c matrixutils.h
        mtwister.c mtwister.h)
