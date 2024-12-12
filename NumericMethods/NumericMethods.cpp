// NumericMethods.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "FredgholmEquation.h"
#include "XLSXWorker.h"
int main()
{
    Fredgholm* fr = new Fredgholm;
    xlsx_worker* excel = new xlsx_worker();
    fr->upper = 3;
    fr->lower = 1;
    fr->step = 10;
    auto sol = fredgholm_solver_by_trapezium(fr);
    excel->write_solution_for_fredgholm_to_excel(L"popa.xls", fr->t, sol);
    auto sol1 = fregholm_solver_by_gauss(fr);
    excel->write_solution_for_fredgholm_to_excel(L"popadr.xls", fr->nodes, sol);
    delete excel;
    delete fr;

}

