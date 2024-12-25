#pragma once
#include "stdafx.h"

class ode_solver {
private:
	xlsx_worker* excel_worker;

	//For stiff ode
	std::vector<double> y_values;
	std::vector<double> z_values;
	std::vector<double> nodes;
	std::vector<double> errors_list;

	//
	std::vector<double> u_values;
	std::vector<double> tau_values;
private:

	double calculate_u_value(double t, double u) {
		return t * std::pow(u, 2);
	}

	double calculate_exact_value(double t) {
		return (-2.0) / (std::pow(t, 2) - 2.0);
	}

	double convert_u_to_implicit(double var_current, double t_current, double step, double step_param) {

		return var_current + (-1.0) * step_param + step * calculate_u_value(t_current, var_current);

	}
	//

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

	double convert_y_to_implicit(double var_current, double z_current, double step, double step_param){

		return var_current + (-1.0) * step_param + step * calculate_y_value(var_current, z_current);
	}

	double convert_z_to_implicit(double var_current, double y_current, double step, double step_param) {

		return var_current + (-1.0) * step_param + step * calculate_z_value(y_current, var_current);
	}



public:
	////
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

		double N = (b - a) / step;

		y_values.resize(N);
		z_values.resize(N);

		y_values[0] = y_s;
		z_values[0] = z_s;

		switch (type) {
			case GEAR_ONE:
				gear_one_solver(step, N);
				calculate_error(L"error_list_gear_one.xls", a, b, step, N);
				break;

			case GEAR_TWO:
				gear_two_solver(step, N);
				calculate_error(L"error_list_gear_two.xls", a, b, step, N);
				break;

			case GEAR_THREE:
				gear_three_solver(step, N);
				calculate_error(L"error_list_gear_three.xls", a, b, step, N);
				break;

			default:
				return;
				
		}
		
		excel_worker->write_solution_for_ode(file_name, N, y_values, z_values);
		
		

		clean_vectors();

	}
	////

	void solve_ode_by_gear(
		int type,
		const wchar_t* file_name,
		double u_s,
		double step,
		double a,
		double b
	)
	{
		double N = (b - a) / step;

		fill_nodes(a, b, step, N);

		u_values.resize(N);

		u_values[0] = 1;

		switch (type) {
		case GEAR_ONE:
			solve_ode_for_u_by_gear_one(step, N);
			calculate_error_for_u(L"error_list_gear_one_u.xls");
			break;

		case GEAR_TWO:
			solve_ode_for_u_by_gear_two(step, N);
			calculate_error_for_u(L"error_list_gear_two_u.xls");
			break;

		case GEAR_THREE:
			solve_ode_for_u_by_gear_three(step, N);
			calculate_error_for_u(L"error_list_gear_three_u.xls");
			break;

		default:
			return;

		}

		excel_worker->write_solution_for_ode(file_name, N, u_values);

		clean_vectors();
	}

	/*
	* @brief Решение ОДУ методом Адомса-Бубубу.
	* @param step - шаг.
	* @param a - начальная граница.
	* @param b - конечная граница.
	*/
	void solve_ode_by_adams_bashfort(
		double step,
		double a,
		double b
		)
	{
		double N = (b - a) / step;

		solve_ode_by_runge_kutta_four(step, N);

		double u_predicator = 0;

		for (size_t i = 3; i < N; i++) {

			u_predicator = u_values[i - 1] + step * (1.5 * calculate_u_value(nodes[i - 1], u_values[i - 1]) - 0.5 * calculate_u_value(nodes[i - 2], u_values[i - 2]));

			u_values[i] = u_values[i - 1] + step * (calculate_u_value(nodes[i], u_predicator) + calculate_u_value(nodes[i - 1], u_values[i - 1])) / 2.0;
		}

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

	void calculate_error_for_u(
		const wchar_t* file_name
		) 
	{
		errors_list.resize(nodes.size());

		for (size_t i = 0; i < nodes.size(); i++) {

			errors_list[i] = std::abs(u_values[i] - calculate_exact_value(nodes[i]));

		}

		excel_worker->write_errors(file_name, nodes, errors_list);

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

	double dihotomy_for_u(
		double a,
		double b,
		double u_current,
		double t_current,
		double step,
		double epsilon = 0.000000001
		)
	{
		double c = 0;
		do
		{
			c = (a + b) / 2;
			if (convert_u_to_implicit(u_current, t_current, step, a) * convert_u_to_implicit(u_current, t_current, step, c) <= 0)
			{
				b = c;
			}
			else
			{
				a = c;
			}
		} while ((b - a) > epsilon);

		return c;
	
	}

	void solve_ode_by_runge_kutta_four(
		double step,
		double N
		)
	{
		fill_nodes(0, 0.1, step, N);

		u_values.resize(N);

		u_values[0] = 1;

		double k1 = 0;
		double k2 = 0;
		double k3 = 0;
		double k4 = 0;

		for (size_t i = 1; i < 3; i++) {

			k1 = calculate_u_value(nodes[i - 1], u_values[i - 1]);
			k2 = calculate_u_value(nodes[i - 1] + 0.5 * step, u_values[i - 1] + 0.5 * step * k1);
			k3 = calculate_u_value(nodes[i - 1] + 0.5 * step, u_values[i - 1] + 0.5 * step * k2);
			k4 = calculate_u_value(nodes[i - 1] + step, u_values[i - 1] + step * k3);

			u_values[i] = u_values[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

		}
	
	}

	void solve_ode_for_u_by_gear_one(
		double step,
		double N
		)
	{
		for (size_t i = 1; i < N; i++) {

			u_values[i] = dihotomy_for_u(-1, 10, u_values[i - 1], nodes[i - 1], step);

		}

	}

	void solve_ode_for_u_by_gear_two(
		double step,
		double N
	)
	{
		solve_ode_for_u_by_gear_one(step, 3);

		for (size_t i = 2; i < N; i++) {

			u_values[i] = dihotomy_for_u(-1, 10, (4.0 / 3.0) * u_values[i - 1] - (1.0 / 3.0) * u_values[i - 2], (2.0 / 3.0) * step, nodes[i-1]);

		}

	}

	void solve_ode_for_u_by_gear_three(
		double step,
		double N
		)
	{
		solve_ode_for_u_by_gear_two(step, N);

		for (size_t i = 3; i < N; i++) {

			u_values[i] = dihotomy_for_u(-1, 10, (18.0 * u_values[i - 1] - 9.0 * u_values[i - 2] + 2 * u_values[i - 3]) / 11.0, 6.0 / 11.0 * step, nodes[i - 1]);

		}
	
	}


	////
	double dihotomy_for_y(
		double a,
		double b,
		double y_current,
		double z_current,
		double step,
		double epsilon = 0.000000001
		)
	{
		double c = 0;
		do
		{
			c = (a + b) / 2;
			if (convert_y_to_implicit(y_current, z_current, step, a) * convert_y_to_implicit(y_current, z_current, step, c) <= 0)
			{
				b = c;
			}
			else
			{
				a = c;
			}
		} while ((b - a) > epsilon);

		return c;
		
	}

	double dihotomy_for_z(
		double a,
		double b,
		double y_current,
		double z_current,
		double step,
		double epsilon = 0.000000001
		)
	{
		double c = 0;
		do
		{
			c = (a + b) / 2;
			if (convert_z_to_implicit(z_current, y_current, step, a) * convert_z_to_implicit(z_current, y_current, step, c) <= 0)
			{
				b = c;
			}
			else
			{
				a = c;
			}
		} while ((b - a) > epsilon);

		return c;

	}

	void gear_one_solver(
		double step,
		double N
		)
	{
		for (size_t i = 1; i < N; i++) {

			y_values[i] = dihotomy_for_y(-1000, 100, y_values[i - 1], z_values[i - 1], step);
			z_values[i] = dihotomy_for_z(-1000, 100, y_values[i - 1], z_values[i - 1], step);

		}

	}

	void gear_two_solver(
		double step,
		double N
		) 
	{
		gear_one_solver(step, 3);

		for (size_t i = 3; i < N; i++) {
			y_values[i] = dihotomy_for_y(-1000, 100000, (4.0 / 3.0) * y_values[i - 1] - (1.0 / 3.0) * y_values[i - 2], z_values[i - 1], (2.0 / 3.0) * step);
			z_values[i] = dihotomy_for_z(-1000, 100000, y_values[i - 1], (4.0 / 3.0) * z_values[i - 1] - (1.0 / 3.0) * z_values[i - 2], (2.0 / 3.0) * step);

		}

	}

	void gear_three_solver(
		double step,
		double N
		) 
	{
		gear_two_solver(step, N);

		for (size_t i = 4; i < N; i++) {
			y_values[i] = dihotomy_for_y(-1000, 100000, (18.0 / 11.0) * y_values[i - 1] - (9.0 / 11.0) * y_values[i - 2] + (2.0 / 11.0) * y_values[i - 3], z_values[i - 1], (6.0 / 11.0) * step);
			z_values[i] = dihotomy_for_z(-1000, 100000, y_values[i - 1], (18.0 / 11.0) * z_values[i - 1] - (9.0 / 11.0) * z_values[i - 2] + (2.0 / 11.0) * y_values[i - 3], (6.0 / 11.0) * step);

		}


	}
	////
private:
	
	void clean_vectors() {
		y_values.clear();
		z_values.clear();
		u_values.clear();
		nodes.clear();
		errors_list.clear();


		y_values.resize(0);
		z_values.resize(0);
		u_values.resize(0);
		nodes.resize(0);
		errors_list.resize(0);
	}

};