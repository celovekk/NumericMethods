#pragma once
#include "stdafx.h"
using namespace libxl;
class xlsx_worker {
public:
	~xlsx_worker() {

		book->release();

	}
	static void write_errors(const wchar_t* file_name, const std::vector<double>& nodes, const std::vector<double>& errors_list) {
		Book* book_writer;

		book_writer = xlCreateBook();
		if (!book_writer)
			return;

		Sheet* sheet = book_writer->addSheet(L"Sheet1");
		if (!sheet)
			return;

		sheet->writeStr(1, 1, L"Nodes");
		sheet->writeStr(1, 2, L"Errors_list");

		for (size_t i = 0; i < nodes.size(); i++) {
			sheet->writeNum(i + 2, 1, nodes[i]);
			sheet->writeNum(i + 2, 2, errors_list[i]);
		}

		book_writer->save(file_name);

		book_writer->release();

	}

	static void write_solution_for_ode(const wchar_t* file_name, double n_count, const std::vector<double>& u_values) {

		Book* book_writer;

		book_writer = xlCreateBook();
		if (!book_writer)
			return;

		Sheet* sheet = book_writer->addSheet(L"Sheet1");
		if (!sheet)
			return;

		sheet->writeStr(1, 1, L"N");
		sheet->writeStr(1, 2, L"u_values");
		

		for (size_t i = 0; i < n_count; i++) {
			sheet->writeNum(i + 2, 1, i);
			sheet->writeNum(i + 2, 2, u_values[i]);
			

		}

		book_writer->save(file_name);

		book_writer->release();
	}

	static void write_solution_for_ode(const wchar_t* file_name, double n_count, const std::vector<double>& y_values, const std::vector<double>& z_values) {

		Book* book_writer;

		book_writer = xlCreateBook();
		if (!book_writer)
			return;

		Sheet* sheet = book_writer->addSheet(L"Sheet1");
		if (!sheet)
			return;

		sheet->writeStr(1, 1, L"N");
		sheet->writeStr(1, 2, L"y_values");
		sheet->writeStr(1, 3, L"z_values");
		
		for (size_t i = 0; i < n_count; i++) {
			sheet->writeNum(i + 2, 1, i);
			sheet->writeNum(i + 2, 2, y_values[i]);
			sheet->writeNum(i + 2, 3, z_values[i]);

		}

		book_writer->save(file_name);

		book_writer->release();
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