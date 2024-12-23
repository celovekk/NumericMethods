// NumericMethods.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "FredgholmEquation.h"
#include "XLSXWorker.h"
#include "ODE.h"
#include "VolterEquation.h"
#include"ODE.h"
int main()
{
    Fredgholm* fr = new Fredgholm;
    Volter* vl = new Volter();
    xlsx_worker* excel = new xlsx_worker();
    ode_solver* od = new ode_solver();
    fr->upper = 3;
    fr->lower = 1;
    fr->step = 10;
    auto sol = fredgholm_solver_by_trapezium(fr);
    excel->write_solution_for_fredgholm_to_excel(L"popa.xls", fr->t, sol);
    auto sol1 = fregholm_solver_by_gauss(fr);
    excel->write_solution_for_fredgholm_to_excel(L"popadr.xls", fr->nodes, sol);

// solve ode
    od->solve_ode_by_implicit_euler(L"implicit_euler.xls", 1, 1, 0.001, 0, 0.1);
    od->solve_ode_by_improved_euler(L"improved_euler.xls", 1, 1, 0.001, 0, 0.1);
    delete od;
    delete excel;
    delete fr;

}

