#pragma once

double LeftRectangleIntegral (double (*F) (double), double a, double b, long parts);
double RightRectangleIntegral (double (*F) (double), double a, double b, long parts);
double MediumRectangleIntegral (double (*F) (double), double a, double b, long parts);
double TrapeziumRectangleIntegral (double (*F) (double), double a, double b, long parts);
double SimpsonRectangleIntegral (double (*F) (double), double a, double b, long parts);
double NewtonIntegral (double (*F) (double), int power, double a, double b, long parts);
