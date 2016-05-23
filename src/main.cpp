/*
 * main.cpp
 *
 *  Created on: May 20, 2016
 *      Author: tobias.neumann
 */

#include "Kernels/SWKernel.h"
#include <iostream>
#include <string>

using std::cout;
using std::string;
using std::endl;

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

	AlignmentKernel * kernel = new SWKernel();

	cout << "Scoring read:\t" << a << endl;
	cout << "Scoring ref:\t" << b << endl;

	float score = kernel->score_alignment(a_s, (float)a.size(), b_s, (float)b.size(),3.0f, 3.0f, 2.0f, 5.0f);

	cout << "Alignment scored " << score << endl;

	delete kernel;
	kernel = 0;

	return 0;
}

