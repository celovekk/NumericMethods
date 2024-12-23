#pragma once
#include "stdafx.h"
class Volter {
	
public:
    std::vector<double> t_values;
    std::vector<double> xt_values;
public:
	double f1(double t, double y1, double y2) {
		return y2; 
	}

	double f2(double t, double y1, double y2) {
		return -4 * y1 + 3 * cos(t);
	}

    double simpson_integral(const std::vector<double>& xVals, const std::vector<double>& tGrid, double t, double tau) {
        int n = t / tau;
        double integral = 0.0;

        for (int i = 1; i < n; i += 2) {
            double s1 = tGrid[i - 1];
            double s2 = tGrid[i];
            double s3 = tGrid[i + 1];
            integral += (s3 - s1) / 6.0 * ((s1 - t) * xVals[i - 1] + 4 * (s2 - t) * xVals[i] + (s3 - t) * xVals[i + 1]);
        }

        return integral;
    }

    void solve_volter(double tau,double N, double T_MAX) {
        std::vector<double> tGrid;
        std::vector<double> xVals;

        // Инициализация сетки t
        for (double t = 0; t <= T_MAX; t += tau) {
            tGrid.push_back(t);
            xVals.push_back(0);  // Начальные условия x(0) = 0
        }

        // Решение интегрального уравнения
        for (size_t i = 1; i < tGrid.size(); ++i) {
            double t = tGrid[i];
            xVals[i] = 4 * simpson_integral(xVals, tGrid, t,tau) + 3 * sin(t);
        }

        t_values = tGrid;
        xt_values = xVals;


    }

	void solve_fde_by_runge_kutta(double t0, double y1_0, double y2_0, double t_end, double h) {
		
       std::vector<double> t_values, x_values;
        double t = t0, y1 = y1_0, y2 = y2_0;


        while (t <= t_end) {
            // Сохраняем текущие значения
            t_values.push_back(t);
            x_values.push_back(y2);  // x(t) = y'(t)

            // Записываем данные в файл
            std::cout << " " << t << " " << y2 << std::endl;

            // Коэффициенты Рунге-Кутты для y1
            double k1_1 = h * f1(t, y1, y2);
            double k2_1 = h * f1(t + h / 2, y1 + k1_1 / 2, y2);
            double k3_1 = h * f1(t + h / 2, y1 + k2_1 / 2, y2);
            double k4_1 = h * f1(t + h, y1 + k3_1, y2);

            // Коэффициенты Рунге-Кутты для y2
            double k1_2 = h * f2(t, y1, y2);
            double k2_2 = h * f2(t + h / 2, y1, y2 + k1_2 / 2);
            double k3_2 = h * f2(t + h / 2, y1, y2 + k2_2 / 2);
            double k4_2 = h * f2(t + h, y1, y2 + k3_2);

            // Обновление y1 и y2
            y1 += (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1) / 6;
            y2 += (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2) / 6;

            // Обновление времени
            t += h;
        }

        
        

	}

    void solve_volter_by() {}

};