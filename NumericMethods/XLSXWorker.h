#pragma once
#include "stdafx.h"
using namespace libxl;
class xlsx_worker {
public:
	~xlsx_worker() {

		book->release();

	}
	void write_solution_for_fredgholm_to_excel(const std::string& file_name, const std::vector<double>& nodes, const std::vector<double>& solution) {
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
				book->save(L"test.xls");
			}
		}
		else {
			return;
		}

		
	}

private:
	Book* book;
};