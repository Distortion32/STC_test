#ifndef HEADER_H
#define HEADER_H

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

#endif // HEADER_H
