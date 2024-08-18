#include <iostream>
#include <cmath>
#include "header.h"
#include <array>

std::array<double, 6> jakobianFunct(const std::array<double, 6> &Coordpoint, const point &d, const point &e, const point &f, const delayStruct &delay)
{
    std::array<double, 9> functArray;
    functArray[0] = sqrt(pow(Coordpoint[0] - d.x,2) + pow(Coordpoint[1] - d.y,2)) - sqrt(pow(Coordpoint[2] - d.x,2) + pow(Coordpoint[3] - d.y,2)) - delay.ADB;
    functArray[1] = sqrt(pow(Coordpoint[0] - d.x,2) + pow(Coordpoint[1] - d.y,2)) - sqrt(pow(Coordpoint[4] - d.x,2) + pow(Coordpoint[5] - d.y,2)) - delay.ADC;
    functArray[2] = sqrt(pow(Coordpoint[2] - d.x,2) + pow(Coordpoint[3] - d.y,2)) - sqrt(pow(Coordpoint[4] - d.x,2) + pow(Coordpoint[5] - d.y,2)) - delay.BDC;

    functArray[3] = sqrt(pow(Coordpoint[0] - e.x,2) + pow(Coordpoint[1] - e.y,2)) - sqrt(pow(Coordpoint[2] - e.x,2) + pow(Coordpoint[3] - e.y,2)) - delay.AEB;
    functArray[4] = sqrt(pow(Coordpoint[0] - e.x,2) + pow(Coordpoint[1] - e.y,2)) - sqrt(pow(Coordpoint[4] - e.x,2) + pow(Coordpoint[5] - e.y,2)) - delay.AEC;
    functArray[5] = sqrt(pow(Coordpoint[2] - e.x,2) + pow(Coordpoint[3] - e.y,2)) - sqrt(pow(Coordpoint[4] - e.x,2) + pow(Coordpoint[5] - e.y,2)) - delay.BEC;

    functArray[6] = sqrt(pow(Coordpoint[0] - f.x,2) + pow(Coordpoint[1] - f.y,2)) - sqrt(pow(Coordpoint[2] - f.x,2) + pow(Coordpoint[3] - f.y,2)) - delay.AFB;
    functArray[7] = sqrt(pow(Coordpoint[0] - f.x,2) + pow(Coordpoint[1] - f.y,2)) - sqrt(pow(Coordpoint[4] - f.x,2) + pow(Coordpoint[5] - f.y,2)) - delay.AFC;
    functArray[8] = sqrt(pow(Coordpoint[2] - f.x,2) + pow(Coordpoint[3] - f.y,2)) - sqrt(pow(Coordpoint[4] - f.x,2) + pow(Coordpoint[5] - f.y,2)) - delay.BFC;
    // подсчитаем якобиан функции F суммы квадратов functArray
    std::array<double, 6> jakobianF;
    jakobianF[0] = 4*((functArray[0] + functArray[1])*(Coordpoint[0] - d.x)/sqrt(pow(Coordpoint[0] - d.x,2) + pow(Coordpoint[1] - d.y,2)) +
                      (functArray[3] + functArray[4])*(Coordpoint[0] - e.x)/sqrt(pow(Coordpoint[0] - e.x,2) + pow(Coordpoint[1] - e.y,2)) +
                      (functArray[6] + functArray[7])*(Coordpoint[0] - f.x)/sqrt(pow(Coordpoint[0] - f.x,2) + pow(Coordpoint[1] - f.y,2)));

    jakobianF[1] = 4*((functArray[0] + functArray[1])*(Coordpoint[1] - d.y)/sqrt(pow(Coordpoint[0] - d.x,2) + pow(Coordpoint[1] - d.y,2)) +
                      (functArray[3] + functArray[4])*(Coordpoint[1] - e.y)/sqrt(pow(Coordpoint[0] - e.x,2) + pow(Coordpoint[1] - e.y,2)) +
                      (functArray[6] + functArray[7])*(Coordpoint[1] - f.y)/sqrt(pow(Coordpoint[0] - f.x,2) + pow(Coordpoint[1] - f.y,2)));

    jakobianF[2] = 4*((functArray[2] - functArray[0])*(Coordpoint[2] - d.x)/sqrt(pow(Coordpoint[2] - d.x,2) + pow(Coordpoint[3] - d.y,2)) +
                      (functArray[5] - functArray[3])*(Coordpoint[2] - e.x)/sqrt(pow(Coordpoint[2] - e.x,2) + pow(Coordpoint[3] - e.y,2)) +
                      (functArray[8] - functArray[6])*(Coordpoint[2] - f.x)/sqrt(pow(Coordpoint[2] - f.x,2) + pow(Coordpoint[3] - f.y,2)));

    jakobianF[3] = 4*((functArray[2] - functArray[0])*(Coordpoint[3] - d.y)/sqrt(pow(Coordpoint[2] - d.x,2) + pow(Coordpoint[3] - d.y,2)) +
                      (functArray[5] - functArray[3])*(Coordpoint[3] - e.y)/sqrt(pow(Coordpoint[2] - e.x,2) + pow(Coordpoint[3] - e.y,2)) +
                      (functArray[8] - functArray[6])*(Coordpoint[3] - f.y)/sqrt(pow(Coordpoint[2] - f.x,2) + pow(Coordpoint[3] - f.y,2)));

    jakobianF[4] = -4*((functArray[1] + functArray[2])*(Coordpoint[4] - d.x)/sqrt(pow(Coordpoint[4] - d.x,2) + pow(Coordpoint[5] - d.y,2)) +
                       (functArray[4] + functArray[5])*(Coordpoint[4] - e.x)/sqrt(pow(Coordpoint[4] - e.x,2) + pow(Coordpoint[5] - e.y,2)) +
                       (functArray[7] + functArray[8])*(Coordpoint[4] - f.x)/sqrt(pow(Coordpoint[4] - f.x,2) + pow(Coordpoint[5] - f.y,2)));

    jakobianF[5] = -4*((functArray[1] + functArray[2])*(Coordpoint[5] - d.y)/sqrt(pow(Coordpoint[4] - d.x,2) + pow(Coordpoint[5] - d.y,2)) +
                       (functArray[4] + functArray[5])*(Coordpoint[5] - e.y)/sqrt(pow(Coordpoint[4] - e.x,2) + pow(Coordpoint[5] - e.y,2)) +
                       (functArray[7] + functArray[8])*(Coordpoint[5] - f.y)/sqrt(pow(Coordpoint[4] - f.x,2) + pow(Coordpoint[5] - f.y,2)));

    // нормируем якобиан
    double normJac = sqrt(pow(jakobianF[0], 2) + pow(jakobianF[1], 2) + pow(jakobianF[2], 2) + pow(jakobianF[3], 2) + pow(jakobianF[4], 2) + pow(jakobianF[5], 2));
    jakobianF[0] = jakobianF[0] / normJac;
    jakobianF[1] = jakobianF[1] / normJac;
    jakobianF[2] = jakobianF[2] / normJac;
    jakobianF[3] = jakobianF[3] / normJac;
    jakobianF[4] = jakobianF[4] / normJac;
    jakobianF[5] = jakobianF[5] / normJac;

    return jakobianF;
}

std::array<double, 6> const_array_multi(const double &LR, const std::array<double, 6> &Jk)
{
    // функция почленного умножения вектора на константу
    std::array<double, 6> lrngRt_JkbnFnct;
    for (int index = 0; index < Jk.size(); ++index)
        lrngRt_JkbnFnct[index] = LR * Jk[index];
    return lrngRt_JkbnFnct;
}

std::array<double, 6> arraySubtraction(const std::array<double, 6> &a, const std::array<double, 6> &b)
{
    // функция вычисления разности векторов
    std::array<double, 6> c;
    for(int index = 0; index < a.size(); ++index)
        c[index] = a[index] - b[index];
    return c;
}

double arrayMod(const std::array<double, 6> &a)
{
    // вычисление модуля вектора
    double mod = 0.0;
    for(int index = 0; index < a.size(); ++index)
        mod += pow(a[index],2);
    return sqrt(mod);
}

double functForCheck(const std::array<double, 6> &Coordpoint, const point &d, const point &e, const point &f, const delayStruct &delay)
{
    std::array<double, 9> functArray;
    functArray[0] = sqrt(pow(Coordpoint[0] - d.x,2) + pow(Coordpoint[1] - d.y,2)) - sqrt(pow(Coordpoint[2] - d.x,2) + pow(Coordpoint[3] - d.y,2)) - delay.ADB;
    functArray[1] = sqrt(pow(Coordpoint[0] - d.x,2) + pow(Coordpoint[1] - d.y,2)) - sqrt(pow(Coordpoint[4] - d.x,2) + pow(Coordpoint[5] - d.y,2)) - delay.ADC;
    functArray[2] = sqrt(pow(Coordpoint[2] - d.x,2) + pow(Coordpoint[3] - d.y,2)) - sqrt(pow(Coordpoint[4] - d.x,2) + pow(Coordpoint[5] - d.y,2)) - delay.BDC;

    functArray[3] = sqrt(pow(Coordpoint[0] - e.x,2) + pow(Coordpoint[1] - e.y,2)) - sqrt(pow(Coordpoint[2] - e.x,2) + pow(Coordpoint[3] - e.y,2)) - delay.AEB;
    functArray[4] = sqrt(pow(Coordpoint[0] - e.x,2) + pow(Coordpoint[1] - e.y,2)) - sqrt(pow(Coordpoint[4] - e.x,2) + pow(Coordpoint[5] - e.y,2)) - delay.AEC;
    functArray[5] = sqrt(pow(Coordpoint[2] - e.x,2) + pow(Coordpoint[3] - e.y,2)) - sqrt(pow(Coordpoint[4] - e.x,2) + pow(Coordpoint[5] - e.y,2)) - delay.BEC;

    functArray[6] = sqrt(pow(Coordpoint[0] - f.x,2) + pow(Coordpoint[1] - f.y,2)) - sqrt(pow(Coordpoint[2] - f.x,2) + pow(Coordpoint[3] - f.y,2)) - delay.AFB;
    functArray[7] = sqrt(pow(Coordpoint[0] - f.x,2) + pow(Coordpoint[1] - f.y,2)) - sqrt(pow(Coordpoint[4] - f.x,2) + pow(Coordpoint[5] - f.y,2)) - delay.AFC;
    functArray[8] = sqrt(pow(Coordpoint[2] - f.x,2) + pow(Coordpoint[3] - f.y,2)) - sqrt(pow(Coordpoint[4] - f.x,2) + pow(Coordpoint[5] - f.y,2)) - delay.BFC;

    double squaredFunc = 0.0;
    for(int index = 0; index < functArray.size(); ++index)
        squaredFunc += pow(functArray[index],2);
    return squaredFunc;
}

int main()
{
    // перечисляем входные данные алгоритма
    // координаты входных точек
    point D{21.0, -19.0};
    point E{-5.0, 16.0};
    point F{2.0, -10.0};
    // координаты искомых точек (они нужны только для подсчёта корректной разности хода сигнала)
    point A{12.0, 7.0};
    point B{10.0, -9.0};
    point C{-14.0, 12.0};

    // разность хода сигнала (можем посчитать точно, так как мы заранее знаем координаты A B C)
    double ADB = sqrt(pow(A.x - D.x,2) + pow(A.y - D.y,2)) - sqrt(pow(B.x - D.x,2) + pow(B.y - D.y,2)); //AD-BD
    double ADC = sqrt(pow(A.x - D.x,2) + pow(A.y - D.y,2)) - sqrt(pow(C.x - D.x,2) + pow(C.y - D.y,2)); //AD-CD
    double BDC = sqrt(pow(B.x - D.x,2) + pow(B.y - D.y,2)) - sqrt(pow(C.x - D.x,2) + pow(C.y - D.y,2)); //BD-CD
    double AEB = sqrt(pow(A.x - E.x,2) + pow(A.y - E.y,2)) - sqrt(pow(B.x - E.x,2) + pow(B.y - E.y,2)); //AE-BE
    double AEC = sqrt(pow(A.x - E.x,2) + pow(A.y - E.y,2)) - sqrt(pow(C.x - E.x,2) + pow(C.y - E.y,2)); //AE-CE
    double BEC = sqrt(pow(B.x - E.x,2) + pow(B.y - E.y,2)) - sqrt(pow(C.x - E.x,2) + pow(C.y - E.y,2)); //BE-CE
    double AFB = sqrt(pow(A.x - F.x,2) + pow(A.y - F.y,2)) - sqrt(pow(B.x - F.x,2) + pow(B.y - F.y,2)); //AF-BF
    double AFC = sqrt(pow(A.x - F.x,2) + pow(A.y - F.y,2)) - sqrt(pow(C.x - F.x,2) + pow(C.y - F.y,2)); //AF-CF
    double BFC = sqrt(pow(B.x - F.x,2) + pow(B.y - F.y,2)) - sqrt(pow(C.x - F.x,2) + pow(C.y - F.y,2)); //BF-CF

    delayStruct delays{ADB, ADC, BDC, AEB, AEC, BEC, AFB, AFC, BFC};

    // задаем точность
    double precision = 0.001;
    // скорость обучения
    double learningRate = 0.2;
    // начальное значение х
    std::array<double, 6> currentPoint{0, 0, 0, 0, 0, 0}; // ax, ay, bx, by, cx, cy
    // значение градиента
    std::array<double, 6> jacob_iteration{0, 0, 0, 0, 0, 0};
    // количество итераций
    int nSteps = 20000;
    int stepIndex = 0;

    // цикл градиентного спуска
    while (stepIndex < nSteps)
    {
        // вычисляем градиент
        jacob_iteration = jakobianFunct(currentPoint, D, E, F, delays);
        // обновление переменного вектора
        std::array<double, 6> newPoint = arraySubtraction(currentPoint, const_array_multi(learningRate, jacob_iteration));
        // проверка условия остановки по сдвигу вектора
        if (arrayMod(arraySubtraction(newPoint, currentPoint)) < precision)
        {
            std::cout << "stepIndex: " << stepIndex << std::endl;
            break;
        }
        // обновление текущего вектора
        currentPoint = newPoint;
        // проверка малости функции минимизации
        if (functForCheck(currentPoint, D, E, F, delays) < precision)
        {
            std::cout << "stepIndex: " << stepIndex << std::endl;
            break;
        }
        stepIndex++;
    }

    std::cout << "A: x = " << currentPoint[0] << "; y = " << currentPoint[1] << std::endl;
    std::cout << "B: x = " << currentPoint[2] << "; y = " << currentPoint[3] << std::endl;
    std::cout << "C: x = " << currentPoint[4] << "; y = " << currentPoint[5] << std::endl;
    std::cout << "Checking the size of the sum of functions squares: " << functForCheck(currentPoint, D, E, F, delays) << std::endl;

    return 0;
}
