#include <iostream>
#include <cmath>
#include "header.h"
#include <array>
//using namespace std;


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
    // Подсчитаем Якобиан функции F суммы квадратов functArray
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

    return jakobianF;
    // ax, ay, bx, by, cx, cy
}



int main()
{
    // Перечисляем входные данные алгоритма
    // Координаты входных точек
    point D{-6.0, 4.0};
    point E{-5.0, 1.0};
    point F{-2.0, -2.0};
    // Координаты искомых точек (они нужны для проверки работы алгоритма и для подсчёта корректной разности хода сигнала)
    point A{5.0, 6.0};
    point B{7.0, 4.0};
    point C{8.0, 1.0};
    // разность хода сигнала (можем посчитать точно, так как мы заранее знаем координаты A B C
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
    double learningRate = 0.9;
    // начальное значение х
    double currentX = 0.9;
    // Значение градиента
    double gradient = 0.0;

    // Цикл градиентного спуска
    while (true)
    {
        gradient = 2 * currentX; //вычисляем градиент
        double newX = currentX - learningRate * gradient; // учитываем learningRate
        // проверяем условие остановки по сдвигу х
        if (fabs(newX - currentX) < precision)
            break;
        currentX = newX; //Обновляем значение х
    }
    std::cout << "Final: " << currentX << std::endl;
    /*
    std::array<double, 6> startPoint{0, 0, 0, 0, 0, 0}; // ax, ay, bx, by, cx, cy
    std::array<double, 6> jacob_iteration;
    jacob_iteration = jakobianFunct(startPoint, D, E, F, delays);

    std::cout << delays.AFB << " " << ADC << " " << BDC << " " << AEB << " " << AEC << " " << BEC << std::endl;
    std::cout << jacob_iteration[0] << " " << jacob_iteration[1] << " " << jacob_iteration[2] << " " << jacob_iteration[3] << " " << jacob_iteration[4] << " " << jacob_iteration[5] << std::endl;
    */
    return 0;
}
