#pragma once
#include "stdafx.h"
class ode_solver{
private:
	std::vector<double> t_values;
	std::vector<double> y_values;
	std::vector<double> z_values;
private:

	double calculate_y_state(double y, double z) {
		return 998 * y + 1998 * z;

	}

	double calculate_z_state(double y, double z) {
		return -999 * y - 1999 * z;
	}
public:
	void solve_ode_by_euler(double step, double a, double b, double y_s, double z_s) {
		double N = (b - a) / step;
		y_values.resize(N);
		z_values.resize(N);

		y_values[0] = y_s;
		z_values[0] = z_s;

		for (size_t i = 1; i < N; i++) {
			y_values[i] = y_values[i - 1] + step * calculate_y_state(y_values[i - 1], z_values[i - 1]);
			z_values[i] = z_values[i - 1] + step * calculate_z_state(y_values[i - 1], z_values[i - 1]);
		}
		excel_worker->write_ode_solution(L"yavniy_euler.xls", N, y_values, z_values);
	}

private:
	xlsx_worker* excel_worker;

};