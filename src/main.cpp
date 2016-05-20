/*
 * main.cpp
 *
 *  Created on: May 20, 2016
 *      Author: tobias.neumann
 */

#include <iostream>
#include <string>

using std::cout;
using std::string;

int main(int argc, char *argv[]) {

	int const RUNS = 1000;

	std::cout << "Startup\n";

	string b = "aaaaaaaaaabbbbbbbaaaaaaaaaa";
	const char * b_s = b.c_str();
	string c = "asdgasdgbmbiabklaklbjasdfja";
	const char * c_s = c.c_str();
	string d = "asdibybkejkbyeubyeubjasdfja";
	const char * d_s = d.c_str();
	string e = "mmmmmmmmmmmmmmmmmmmmjasdfja";
	const char * e_s = e.c_str();

	string a = "aaaaaaaaaaaaaaaaaaaa";
	const char *a_s = a.c_str();
	string a1 = "asdgasdgbmaaaaaaaaaa";
	const char *a1_s = a1.c_str();
	string a2 = "aavdvgasejkbyeubyaaa";
	const char *a2_s = a2.c_str();
	string a3 = "asdgasdgbmaaaaaaaaaa";
	const char *a3_s = a3.c_str();
	std::cout << "Hello world" << std::endl;
	return 0;
}

