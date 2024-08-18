#ifndef HEADER_H
#define HEADER_H
#include <array>

// Структура для хранения координат точки
struct point{
    double x;
    double y;
};

// Структура для хранения разности хода сигнала
struct delayStruct{
    double ADB; //AD-BD
    double ADC; //AD-CD
    double BDC; //BD-CD
    double AEB; //AE-BE
    double AEC; //AE-CE
    double BEC; //BE-CE
    double AFB; //AF-BF
    double AFC; //AF-CF
    double BFC; //BF-CF
};

std::array<double, 6> jakobianFunct(const std::array<double, 6> &Coordpoint, const point &d, const point &e, const point &f, const delayStruct &delay);
// функция подсчета вектора частных производных
std::array<double, 6> const_array_multi(const double &LR, const std::array<double, 6> &Jk);
// функция почленного умножения вектора на константу
std::array<double, 6> arraySubtraction(const std::array<double, 6> &a, const std::array<double, 6> &b);
// функция вычисления разности векторов
double arrayMod(const std::array<double, 6> &a);
// вычисление модуля вектора
double functForCheck(const std::array<double, 6> &Coordpoint, const point &d, const point &e, const point &f, const delayStruct &delay);
// метод для проверки малости функции, которую нужно минимизировать

#endif // HEADER_H
