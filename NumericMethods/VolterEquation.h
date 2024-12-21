#pragma once
#include "stdafx.h"
class Volter {
	
public:
	double f1(double t, double y1, double y2) {
		return y2; 
	}

	double f2(double t, double y1, double y2) {
		return -4 * y1 + 3 * cos(t);
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

};