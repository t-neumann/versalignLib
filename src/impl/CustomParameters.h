#include "AlignmentParameters.h"

#include <string.h>
#include <string>

class CustomParameters : public AlignmentParameters {

	public :
		virtual int param_int(char const * const key) {
			if (strcmp(key, "score_match") == 0) {
				return score_match;
			} else if (strcmp(key, "score_mismatch") == 0) {
				return score_mismatch;
			} else if (strcmp(key, "score_gap_read") == 0) {
				return score_gap_read;
			} else if (strcmp(key, "score_gap_ref") == 0) {
				return score_gap_ref;
			} else if (strcmp(key, "read_length") == 0) {
				return read_length;
			} else if (strcmp(key, "ref_length") == 0) {
				return ref_length;
			}
			std::string exception("Unknown int parameter: ");
			exception += key;
			throw exception.c_str();
		}

		virtual bool has_key(char const * const key) {
			if (strcmp(key, "score_match") == 0) {
				return true;
			} else if (strcmp(key, "score_mismatch") == 0) {
				return true;
			} else if (strcmp(key, "score_gap_read") == 0) {
				return true;
			} else if (strcmp(key, "score_gap_ref") == 0) {
				return true;
			} else if (strcmp(key, "read_length") == 0) {
				return true;
			} else if (strcmp(key, "ref_length") == 0) {
				return true;
			}
			return false;
		}

		int read_length = 0;
		int ref_length = 0;

	private :

		short score_match = 2;
		short score_mismatch = -1;
		short score_gap_read = -3;
		short score_gap_ref = -3;
};
