#pragma once
#include "stdafx.h"
using namespace libxl;
class xlsx_worker {
public:
	~xlsx_worker() {

		book->release();

	}

	static void write_ode_solution(const wchar_t* file_name, double N, const std::vector<double>& y_values, const std::vector<double>& z_values) {
		Book* booka = nullptr;
		if (booka) 
			booka->release();

		booka = xlCreateBook();

		if (booka) {
			Sheet* sheet = booka->addSheet(L"Sheet1");
			if (sheet) {
				sheet->writeStr(1, 1, L"steps");
				sheet->writeStr(1, 2, L"y_val");
				sheet->writeStr(1, 3, L"z_val");

				for (size_t i = 0; i < N; i++) {
					sheet->writeNum(i + 2, 1, i);
					sheet->writeNum(i + 2, 2, y_values[i]);
					sheet->writeNum(i + 2, 3, z_values[i]);
				}
				booka->save(file_name);
			}
			else { return; }
		}
	}

	void write_solution_for_fredgholm_to_excel(const wchar_t* file_name, const std::vector<double>& nodes, const std::vector<double>& solution) {
		if (book) {
			book->release();
		}
		book = xlCreateBook();

		if (book) {
			Sheet* sheet = book->addSheet(L"Sheet1");
			if (sheet) {
				sheet->writeStr(1, 1, L"Nodes");
				sheet->writeStr(1, 2, L"Solution");
				
				for (size_t i = 0; i < nodes.size(); i++) {
					sheet->writeNum(i + 2, 1, nodes[i]);
					sheet->writeNum(i + 2, 2, solution[i]);
				}
				book->save(file_name);
			}
		}
		else {
			return;
		}

		
	}

private:
	Book* book;
};