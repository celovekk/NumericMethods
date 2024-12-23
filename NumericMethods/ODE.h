#pragma once
#include "stdafx.h"

class ode_solver {
private:
	xlsx_worker* excel_worker;

	std::vector<double> y_values;
	std::vector<double> z_values;

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

	double analytical_y(double z) {
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

		for(size_t i = 0; i < N; i++){

			y_values[i] = y_values[i - 1] + step * calculate_y_value(y_values[i - 1], z_values[i - 1]);
			z_values[i] = z_values[i - 1] + step * calculate_z_value(y_values[i - 1], z_values[i - 1]);

		}

		excel_worker->write_solution_for_ode(file_name, N, y_values, z_values);

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
	* @param b - êîíå÷íîå конечная граница.
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
	* @param b - êîíå÷íîå конечная граница.
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

		excel_worker->write_ode_solution(L"yavniy_euler.xls", N, y_values, z_values);
	}

private:
	
	void calculate_error(
		std::function<double(double)> analytical_val_y,
		std::function<double(double)> analytical_val_z,

		) {}

private:
	void clean_vectors() {
		y_values.clear();
		z_values.clear();

		y_values.resize(0);
		z_values.resize(0);
	}

	

};