// NumericMethods.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "FredgholmEquation.h"
int main()
{
    Fredgholm fr ;
    fr.upper = 3;
    fr.lower = 1;
    fr.step = 10;
    auto ten = fregholm_solver_by_gauss(&fr);
    for (double sol : ten) {
        std::cout << sol << " ";
    }
    for (int i = 0; i < ten.size(); ++i) {
        double t = fr.lower + i * (fr.upper - fr.lower) / (ten.size() - 1);
        std::cout << "x(" << t << ") = " << ten[i] << "\n";
    }

}

