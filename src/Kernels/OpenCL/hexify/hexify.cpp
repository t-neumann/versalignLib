#include <iostream>
#include <fstream>
#include <string>

using std::ifstream;
using std::endl;
using std::string;
using std::cout;
using std::ofstream;
using std::hex;

//#include <unistd.h>

#if defined(WIN32) || defined(_WIN32)
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif

int main(int argc, char *argv[]) {

//	char cwd[1024];
//	if (getcwd(cwd, sizeof(cwd)) != NULL)
//		fprintf(stdout, "Current working dir: %s\n", cwd);
//	else
//		perror("getcwd() error");

	string file (argv[argc - 1]);

	string sourcename = file;

	const size_t last_slash_idx = sourcename.find_last_of(PATH_SEPARATOR);
	if (std::string::npos != last_slash_idx)
	{
		sourcename.erase(0, last_slash_idx + 1);
	}

//	cout << "OpenCL file: " << file << endl;

	ifstream clfile ((file + ".cl").c_str());

	if (clfile.is_open()) {
		ofstream hfile ((file + ".h").c_str());
		if (hfile.is_open()) {
			string line;

			hfile << "char const " << sourcename << "[] = {\n";
			while(getline(clfile, line)) {
				for (int i = 0; i < line.length(); ++i) {
					hfile << "0x" << std::hex << (int)line[i] << ", ";
				}
				hfile << "0x" << std::hex << (int)'\n' << ", " ;
			}
			hfile << "0x" << std::hex << (int)'\0' << "};\n";
			hfile.close();
		} else {
			cout << "Unable to open output file" << endl;
		}
		clfile.close();
	} else {
		cout << "Unable to open cl file" << endl;
	}
	return 0;
}

