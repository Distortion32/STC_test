#include <iostream>
#include <cmath>
#include "header.h"
using namespace std;

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

    cout << ADB << " " << ADC << " " << BDC << " " << AEB << " " << AEC << " " << BEC << endl;
    return 0;
}
