// NumericMethods.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "FredgholmEquation.h"
#include "XLSXWorker.h"
#include "ODE.h"
#include "VolterEquation.h"
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
    vl->solve_volter(0.25, 20, 5.0);
    excel->write_solution_for_fredgholm_to_excel(L"volter.xls", vl->t_values, vl->xt_values);
    od->solve_ode_by_euler(0.001, 0, 0.1, 1, 1);
    delete od;
    delete vl;
    delete excel;
    delete fr;

}

