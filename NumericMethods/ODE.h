#pragma once
#include "stdafx.h"

class ode_solver {
private:
	xlsx_worker* excel_worker;

	std::vector<double> y_values;
	std::vector<double> z_values;
	std::vector<double> nodes;
	std::vector<double> errors_list;

private:

	double calculate_y_value(double y, double z){
		return 998 * y + 1998 * z;
	}

	double calculate_z_value(double y, double z) {
		return -999 * y - 1999 * z;
	}

	double analytical_y(double y) {
		return exp(-1000 * y) * (-3 + 4 * exp(999 * y));
	}

	double analytical_z(double z) {
		return exp(-1000 * z) * (3 - 2 * exp(999 * z));
	}

public:

	/*
	* @brief Решение ОДУ явным методом Эйлера.
	*
	* @param file_name - путь, куда сохранить ecxel файл.
	* @param y_s - Y(0).
	* @param z_s - Z(0).
	* @param step - шаг.
	* @param a - начальная граница.
	* @param b - êîíå÷íîå конечная граница.
	*/
	void solve_ode_by_euler(
		const wchar_t* file_name,
		double y_s,
		double z_s,
		double step,
		double a,
		double b
	) 
	{
		double N = (b - a) / step;

		y_values.resize(N);
		z_values.resize(N);

		y_values[0] = y_s;
		z_values[0] = z_s;

		for(size_t i = 1; i < N; i++){

			y_values[i] = y_values[i - 1] + step * calculate_y_value(y_values[i - 1], z_values[i - 1]);
			z_values[i] = z_values[i - 1] + step * calculate_z_value(y_values[i - 1], z_values[i - 1]);

		}

		excel_worker->write_solution_for_ode(file_name, N, y_values, z_values);

		calculate_error(L"errors_list_for_euler.xls", a, b, step, N);

		clean_vectors();
	}

	/*
	* @brief Решение ОДУ неявным методом Эйлера.
	*
	* @param file_name - путь, куда сохранить ecxel файл.
	* @param y_s - Y(0).
	* @param z_s - Z(0).
	* @param step - шаг.
	* @param a - начальная граница.
	* @param b - конечная граница.
	*/
	void solve_ode_by_implicit_euler(
		const wchar_t* file_name,
		double y_s,
		double z_s,
		double step,
		double a,
		double b
	)
	{
		double N = (b - a) / step;

		y_values.resize(N);
		z_values.resize(N);

		y_values[0] = y_s;
		z_values[0] = z_s;

		for (size_t i = 1; i < N; i++) {
			y_s = y_values[i - 1] + step * calculate_y_value(y_values[i - 1], z_values[i - 1]);
			z_s = z_values[i - 1] + step * calculate_z_value(y_values[i - 1], z_values[i - 1]);

			y_values[i] = y_values[i - 1] + step * (calculate_y_value(y_values[i - 1], z_values[i - 1]) + calculate_y_value(y_s, z_s)) / 2;
			z_values[i] = z_values[i - 1] + step * (calculate_z_value(y_values[i - 1], z_values[i - 1]) + calculate_z_value(y_s, z_s)) / 2;
		}

		excel_worker->write_solution_for_ode(file_name, N, y_values, z_values);

		clean_vectors();

	}

	/*
	* @brief Решение ОДУ улучшенным методом Эйлера второго порядка.
	*
	* @param file_name - путь, куда сохранить ecxel файл.
	* @param y_s - Y(0).
	* @param z_s - Z(0).
	* @param step - шаг.
	* @param a - начальная граница.
	* @param b - конечная граница.
	*/
	void solve_ode_by_improved_euler(
		const wchar_t* file_name,
		double y_s,
		double z_s,
		double step,
		double a,
		double b
		)
	{
		double N = (b - a) / step;

		y_values.resize(N);
		z_values.resize(N);

		y_values[0] = y_s;
		z_values[0] = z_s;

		for (size_t i = 1; i < N; i++) {
			y_values[i] = y_values[i - 1] + step * calculate_y_value(
				y_values[i - 1] + step / 2 * calculate_y_value(y_values[i - 1], z_values[i - 1]),
				z_values[i - 1] + step / 2 * calculate_z_value(y_values[i - 1], z_values[i - 1]));

			z_values[i] = z_values[i - 1] + step * calculate_z_value(
				y_values[i - 1] + step / 2 * calculate_y_value(y_values[i - 1], z_values[i - 1]),
				z_values[i - 1] + step / 2 * calculate_z_value(y_values[i - 1], z_values[i - 1]));
		}

		excel_worker->write_solution_for_ode(file_name, N, y_values, z_values);

		clean_vectors();

	}

	/*
	* @brief Решение ОДУ улучшенным методами Гира 1, 2, 4 порядка.
	*
	* @param type - 
	* GEAR_ONE: метод Гира 1 порядка,
	* GEAR_TWO: метод Гира 2 порядка,
	* GEAR_FOUR: метод Гира 4 порядка.
	* @param file_name - путь, куда сохранить ecxel файл.
	* @param y_s - Y(0).
	* @param z_s - Z(0).
	* @param step - шаг.
	* @param a - начальная граница.
	* @param b - конечная граница.
	*/
	void solve_ode_by_gear(
		int type,
		const wchar_t* file_name,
		double y_s,
		double z_s,
		double step,
		double a,
		double b
		) 
	{
		int iter = 0;

		double y_corr, z_corr;	
		double N = (b - a) / step;

		y_values.resize(N);
		z_values.resize(N);

		y_values[0] = y_s;
		z_values[0] = z_s;

		switch (type) {
			case GEAR_ONE:
				for (size_t i = 1; i < N; i++) {

					y_values[i] = y_values[i - 1] + step * calculate_y_value(y_values[i - 1], z_values[i - 1]);
					z_values[i] = z_values[i - 1] + step * calculate_z_value(y_values[i - 1], z_values[i - 1]);

				}
				break;

			case GEAR_TWO:
				for (size_t i = 1; i < N; i++) {					
					if (i != 1) {
						y_s = y_values[i - 1] + step * calculate_y_value(y_values[i - 1], z_values[i - 1]);
						z_s = z_values[i - 1] + step * calculate_z_value(y_values[i - 1], z_values[i - 1]);

						y_values[i] = (4.0 * y_values[i - 1] - y_values[i - 2] + 2 * step * calculate_y_value(y_s, z_s)) / 3;
						z_values[i] = (4.0 * z_values[i - 1] - z_values[i - 2] + 2 * step * calculate_z_value(y_s, z_s)) / 3;

					}
					else {
						y_values[i] = y_values[i - 1] + step * calculate_y_value(y_values[i - 1], z_values[i - 1]);
						z_values[i] = z_values[i - 1] + step * calculate_z_value(y_values[i - 1], z_values[i - 1]);
					}
					
					
				}
				break;

			case GEAR_FOUR:
				for (size_t i = 1; i < 4; i++) {
					y_values[i] = y_values[i - 1] + step * calculate_y_value(y_values[i - 1], z_values[i - 1]);
					z_values[i] = z_values[i - 1] + step * calculate_z_value(y_values[i - 1], z_values[i - 1]);

				}

				for (size_t i = 4; i < N; i++) {
					y_s = y_values[i - 1] + step * calculate_y_value(y_values[i - 1], z_values[i - 1]);
					z_s = z_values[i - 1] + step * calculate_z_value(y_values[i - 1], z_values[i - 1]);

					y_values[i] = (48 * y_values[i - 1] - 36 * y_values[i - 2] + 16 * y_values[i - 3] + 3 * y_values[i - 4] + 12 * step * calculate_y_value(y_s, z_s)) / 25;
					z_values[i] = (48 * z_values[i - 1] - 36 * z_values[i - 2] + 16 * z_values[i - 3] + 3 * z_values[i - 4] + 12 * step * calculate_z_value(y_s, z_s)) / 25;

				}
				break;

			default:
				return;
				

		}
		
		excel_worker->write_solution_for_ode(file_name, N, y_values, z_values);

		clean_vectors();

	}

private:
	std::vector<double> fill_nodes(double a, double b, double h, double n) {
		nodes.resize(n);

		for (size_t i = 0; i < n; i++) {
			a = i * h;
			nodes[i] = a;
		}

		return nodes;
	}

	void calculate_error(
		const wchar_t* file_name,
		double a,
		double b,
		double h,
		double n
		) 
	{
		nodes = fill_nodes(a, b, h, n);

		errors_list.resize(nodes.size());

		double def_y = 0;
		double def_z = 0;

		for (size_t i = 0; i < nodes.size(); i++) {

			def_y = abs(analytical_y(nodes[i]) - y_values[i]);
			def_z = abs(analytical_z(nodes[i]) - z_values[i]);
			
			errors_list[i] = std::max(def_y, def_z);
		}

		excel_worker->write_errors(file_name, nodes, errors_list);	
	
	}

private:
	
	void clean_vectors() {
		y_values.clear();
		z_values.clear();
		nodes.clear();
		errors_list.clear();

		y_values.resize(0);
		z_values.resize(0);
		nodes.resize(0);
		errors_list.resize(0);
	}

	

};