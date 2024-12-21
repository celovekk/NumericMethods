// NumericMethods.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "FredgholmEquation.h"
#include "XLSXWorker.h"
#include "VolterEquation.h"
int main()
{
    Fredgholm* fr = new Fredgholm;
    Volter* vl = new Volter();
    xlsx_worker* excel = new xlsx_worker();
    fr->upper = 3;
    fr->lower = 1;
    fr->step = 10;
    auto sol = fredgholm_solver_by_trapezium(fr);
    //excel->write_solution_for_fredgholm_to_excel("popa.xls", fr->weights, sol);
    vl->solve_fde_by_runge_kutta(0, 0, 0, 5.0, 0.25);
    delete vl;
    delete excel;
    delete fr;

}

