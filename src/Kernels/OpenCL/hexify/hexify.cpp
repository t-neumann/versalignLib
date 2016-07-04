#include <iostream>
#include <fstream>
#include <string>

using std::ifstream;
using std::endl;
using std::string;
using std::cout;
using std::ofstream;
using std::hex;

int main(int argc, char *argv[]) {

	string name (argv[argc - 1]);

	ifstream clfile ((name + ".cl").c_str());

	if (clfile.is_open()) {
		ofstream hfile ((name + ".h").c_str());
		if (hfile.is_open()) {
			string line;

			hfile << "char const " << name << "[] = {\n";
			while(getline(clfile, line)) {
				for (int i = 0; i < line.length(); ++i) {
					hfile << "0x" << std::hex << (int)line[i] << ", ";
				}
				hfile << "0x" << std::hex << (int)'\n' << ", " ;
				}
			hfile << "0x" << std::hex << (int)'\0' << "};\n";
		} else {
			cout << "Unable to open output file" << endl;
		}
	} else {
		cout << "Unable to open cl file" << endl;
	}
	return 0;
}

