__kernel void score_alignment_smith_waterman(__global const char* read, __global const char* ref, __global short * results) {

	// Offset read and ref pointers
	read = read + g_id * read_length * VSIZE;
	ref = ref + g_id * ref_length * VSIZE;

	short16 max_score = v_null;

	// Initialize matrix
	short16 matrix [(ref_length + 1) * SCORING_ROWS];

	// offset for first and second row
	int prev_row = 0;
	int cur_row = 1;

	for (int ref_pos = 0; ref_pos < ref_length + 1; ++ref_pos) {
		matrix[ref_pos] = v_null;
	}

	for (int read_pos = 0; read_pos < read_length; ++read_pos) {
	
		char16 read_cache = (char16)(
			char_to_score[read[read_pos]],
			char_to_score[read[read_pos + read_length]],
			char_to_score[read[read_pos + 2 * (read_length)]],
			char_to_score[read[read_pos + 3 * (read_length)]],
			char_to_score[read[read_pos + 4 * (read_length)]],
			char_to_score[read[read_pos + 5 * (read_length)]],
			char_to_score[read[read_pos + 6 * (read_length)]],
			char_to_score[read[read_pos + 7 * (read_length)]],
			char_to_score[read[read_pos + 8 * (read_length)]],
			char_to_score[read[read_pos + 9 * (read_length)]],
			char_to_score[read[read_pos + 10 * (read_length)]],
			char_to_score[read[read_pos + 11 * (read_length)]],
			char_to_score[read[read_pos + 12 * (read_length)]],
			char_to_score[read[read_pos + 13 * (read_length)]],
			char_to_score[read[read_pos + 14 * (read_length)]],
			char_to_score[read[read_pos + 15 * (read_length)]]
		);

		matrix[cur_row * (ref_length + 1)] = v_null;

		//printf("Read base =%c and score=%i\n", read[read_pos],read_cache.s0);
		//printf("Read2 base =%c and score2=%i\n", read[read_pos+ read_length],read_cache.s1);

		for (int ref_pos = 0; ref_pos < ref_length; ++ref_pos) {
		
			short16 base_cmp = (short16)(base_score[read_cache.s0][char_to_score[ref[ref_pos]]],
										 base_score[read_cache.s1][char_to_score[ref[ref_pos + ref_length]]],
										 base_score[read_cache.s2][char_to_score[ref[ref_pos + 2 *ref_length]]],
										 base_score[read_cache.s3][char_to_score[ref[ref_pos + 3 *ref_length]]],
										 base_score[read_cache.s4][char_to_score[ref[ref_pos + 4 *ref_length]]],
										 base_score[read_cache.s5][char_to_score[ref[ref_pos + 5 *ref_length]]],
										 base_score[read_cache.s6][char_to_score[ref[ref_pos + 6 *ref_length]]],
										 base_score[read_cache.s7][char_to_score[ref[ref_pos + 7 *ref_length]]],
										 base_score[read_cache.s8][char_to_score[ref[ref_pos + 8 *ref_length]]],
										 base_score[read_cache.s9][char_to_score[ref[ref_pos + 9 *ref_length]]],
										 base_score[read_cache.sA][char_to_score[ref[ref_pos + 10 *ref_length]]],
										 base_score[read_cache.sB][char_to_score[ref[ref_pos + 11 *ref_length]]],
										 base_score[read_cache.sC][char_to_score[ref[ref_pos + 12 *ref_length]]],
										 base_score[read_cache.sD][char_to_score[ref[ref_pos + 13 *ref_length]]],
										 base_score[read_cache.sE][char_to_score[ref[ref_pos + 14 *ref_length]]],
										 base_score[read_cache.sF][char_to_score[ref[ref_pos + 15 *ref_length]]]
			);

			// up = matrix[prev_row * (ref_length + 1) + ref_pos + 1];
			// diag = matrix[prev_row * (ref_length + 1) + ref_pos];
			// left = matrix[current_row * (ref_length + 1) + ref_pos];
			
			//short16 up = matrix[prev_row * (ref_length + 1) + ref_pos + 1] + v_score_gap_ref;
			//short16 left = matrix[cur_row * (ref_length + 1) + ref_pos] + v_score_gap_read;
			//short16 diag = matrix[prev_row * (ref_length + 1) + ref_pos] + base_cmp;

			short16 cur = max(matrix[prev_row * (ref_length + 1) + ref_pos + 1] + v_score_gap_ref,
						  max(matrix[cur_row * (ref_length + 1) + ref_pos] + v_score_gap_read,
						  max(matrix[prev_row * (ref_length + 1) + ref_pos] + base_cmp,
						  v_null)));

			//printf("UP =%hi ", up.s0);
			//printf("LEFT =%hi ", left.s0);
			//printf("DIAG =%hi ", diag.s0);
			//printf("CUR =%hi\n", cur.s0);

			matrix[cur_row * (ref_length + 1) + ref_pos + 1] =  cur;

			max_score = max(max_score, cur);

		}
		//printf("\nMAX_SCORE = %hi\n",max_score.s0);
		prev_row = cur_row++;
		cur_row %= SCORING_ROWS;
	}
	vstore16(max_score, g_id, results);
}

__kernel void score_alignment_needleman_wunsch(__global const char* read, __global const char* ref, __global short * results) {

	// Offset read and ref pointers
	read = read + g_id * read_length * VSIZE;
	ref = ref + g_id * ref_length * VSIZE;

	short16 max_score = (short16)(SHORT_MIN);

	// Initialize matrix
	short16 matrix [(ref_length + 1) * SCORING_ROWS];

	// offset for first and second row
	int prev_row = 0;
	int cur_row = 1;

	for (int ref_pos = 0; ref_pos < ref_length + 1; ++ref_pos) {
		matrix[ref_pos] = v_null;
	}

	for (int read_pos = 0; read_pos < read_length; ++read_pos) {

		char16 read_cache = (char16)(
				char_to_score[read[read_pos]],
				char_to_score[read[read_pos + read_length]],
				char_to_score[read[read_pos + 2 * (read_length)]],
				char_to_score[read[read_pos + 3 * (read_length)]],
				char_to_score[read[read_pos + 4 * (read_length)]],
				char_to_score[read[read_pos + 5 * (read_length)]],
				char_to_score[read[read_pos + 6 * (read_length)]],
				char_to_score[read[read_pos + 7 * (read_length)]],
				char_to_score[read[read_pos + 8 * (read_length)]],
				char_to_score[read[read_pos + 9 * (read_length)]],
				char_to_score[read[read_pos + 10 * (read_length)]],
				char_to_score[read[read_pos + 11 * (read_length)]],
				char_to_score[read[read_pos + 12 * (read_length)]],
				char_to_score[read[read_pos + 13 * (read_length)]],
				char_to_score[read[read_pos + 14 * (read_length)]],
				char_to_score[read[read_pos + 15 * (read_length)]]
		);

		matrix[cur_row * (ref_length + 1)] = v_null;

		//printf("Read base =%c and score=%i\n", read[read_pos],read_cache.s0);
		//printf("Read2 base =%c and score2=%i\n", read[read_pos+ read_length],read_cache.s1);

		for (int ref_pos = 0; ref_pos < ref_length; ++ref_pos) {

			short16 base_cmp = (short16)(base_score[read_cache.s0][char_to_score[ref[ref_pos]]],
					base_score[read_cache.s1][char_to_score[ref[ref_pos + ref_length]]],
					base_score[read_cache.s2][char_to_score[ref[ref_pos + 2 *ref_length]]],
					base_score[read_cache.s3][char_to_score[ref[ref_pos + 3 *ref_length]]],
					base_score[read_cache.s4][char_to_score[ref[ref_pos + 4 *ref_length]]],
					base_score[read_cache.s5][char_to_score[ref[ref_pos + 5 *ref_length]]],
					base_score[read_cache.s6][char_to_score[ref[ref_pos + 6 *ref_length]]],
					base_score[read_cache.s7][char_to_score[ref[ref_pos + 7 *ref_length]]],
					base_score[read_cache.s8][char_to_score[ref[ref_pos + 8 *ref_length]]],
					base_score[read_cache.s9][char_to_score[ref[ref_pos + 9 *ref_length]]],
					base_score[read_cache.sA][char_to_score[ref[ref_pos + 10 *ref_length]]],
					base_score[read_cache.sB][char_to_score[ref[ref_pos + 11 *ref_length]]],
					base_score[read_cache.sC][char_to_score[ref[ref_pos + 12 *ref_length]]],
					base_score[read_cache.sD][char_to_score[ref[ref_pos + 13 *ref_length]]],
					base_score[read_cache.sE][char_to_score[ref[ref_pos + 14 *ref_length]]],
					base_score[read_cache.sF][char_to_score[ref[ref_pos + 15 *ref_length]]]
			);

			// up = matrix[prev_row * (ref_length + 1) + ref_pos + 1];
			// diag = matrix[prev_row * (ref_length + 1) + ref_pos];
			// left = matrix[current_row * (ref_length + 1) + ref_pos];

			//short16 up = matrix[prev_row * (ref_length + 1) + ref_pos + 1] + v_score_gap_ref;
			//short16 left = matrix[cur_row * (ref_length + 1) + ref_pos] + v_score_gap_read;
			//short16 diag = matrix[prev_row * (ref_length + 1) + ref_pos] + base_cmp;

			short16 cur = max(matrix[prev_row * (ref_length + 1) + ref_pos + 1] + v_score_gap_ref,
					max(matrix[cur_row * (ref_length + 1) + ref_pos] + v_score_gap_read,
							matrix[prev_row * (ref_length + 1) + ref_pos] + base_cmp));

			//printf("UP =%hi ", up.s0);
			//printf("LEFT =%hi ", left.s0);
			//printf("DIAG =%hi ", diag.s0);


			matrix[cur_row * (ref_length + 1) + ref_pos + 1] =  cur;

		}

		max_score = max(max_score, matrix[cur_row * (ref_length + 1) + ref_length]);

		//printf("\nMAX_SCORE = %hi\n",max_score.s0);
		prev_row = cur_row++;
		cur_row %= SCORING_ROWS;
	}

	for (int ref_pos = 0; ref_pos < ref_length + 1; ++ref_pos) {
		max_score = max(max_score, matrix[prev_row * (ref_length + 1) + ref_pos]);
	}

	vstore16(max_score, g_id, results);
}
