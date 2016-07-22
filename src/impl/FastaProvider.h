/*
 * FastaProvider.h
 *
 *  Created on: Jul 22, 2016
 *      Author: tobias.neumann
 */

#ifndef FASTAPROVIDER_H
#define FASTAPROVIDER_H

#include <fstream>
#include <string>
#include <iostream>
#include <vector>

using std::ifstream;
using std::cout;
using std::getline;
using std::string;
using std::vector;

class FastaProvider {

public:
	FastaProvider() {
	}

	vector<const char *> parse_fasta(string const & filename) {

		vector<const char*> values;

		ifstream in(filename.c_str());

		if (!in.good()) {
			std::cerr << "Error opening " + filename << ".\n";
		} else {
			std::string line, name, content;
			while (getline(in, line).good()) {
				if (line.empty() || line[0] == '>') {
					if (!name.empty()) {
						//char * charSeq = new char[content.length() + 1];
						//strncpy(charSeq, content.c_str(), content.length() + 1);
						//values.push_back(charSeq);
						values.push_back(strdup(content.c_str()));
						name.clear();
					}
					if (!line.empty()) {
						name = line.substr(1);
					}
					content.clear();
				} else if (!name.empty()) {
					if (line.find(' ') != std::string::npos) {
						name.clear();
						content.clear();
					} else {
						content += line;
					}
				}
			}
			if (!name.empty()) {
				//char * charSeq = new char[content.length() + 1];
				//strncpy(charSeq, content.c_str(), content.length() + 1);
				//values.push_back(charSeq);
				values.push_back(strdup(content.c_str()));
			}
		}
		return values;
	}



	~FastaProvider() {
	}

};

#endif /* FASTAPROVIDER_H */
