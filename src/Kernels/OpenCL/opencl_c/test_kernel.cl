#pragma OPENCL EXTENSION cl_amd_printf : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store: enable

#ifdef __NVIDIA__
#define __ocl_constant __constant const
#else
#define __ocl_constant __constant
#endif

// OCL indices

#define g_id get_global_id(0)
#define l_id get_local_id(0)

// matrix pointers
#define UP 1
#define LEFT 2
#define DIAG 3
#define START 0

// base translation
#define ASCII_ALPHABET 256
#define SCORE_CASE 6

__ocl_constant char char_to_score[ASCII_ALPHABET] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 3, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 5, 0,
    0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 3, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 5, 0,
    0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// SSE = 8 short per instruction
// AVX = 16 short per instruction
// (trying 16-vectorized variables now)
#define VSIZE 16

#define SHORT_MIN -32767

#define SCORING_ROWS 2

// Scores

__ocl_constant short16 v_null = (short16)(0);
__ocl_constant short16 v_start = (short16)(START);
__ocl_constant short16 v_up = (short16)(UP);
__ocl_constant short16 v_left = (short16)(LEFT);
__ocl_constant short16 v_diag = (short16)(DIAG);

__ocl_constant char16 v_nullchar = (char16)('\0');

__ocl_constant short16 v_score_gap_read = (short16)(score_gap_read);
__ocl_constant short16 v_score_gap_ref = (short16)(score_gap_ref);
__ocl_constant short16 v_score_match = (short16)(score_match);
__ocl_constant short16 v_score_mismatch = (short16)(score_mismatch);


__ocl_constant short base_score[SCORE_CASE][SCORE_CASE]= {
    // non ATGCN
    {0,0,0,0,0,0},
    // A
    {0,score_match,score_mismatch,score_mismatch,score_mismatch,0},
    // T
    {0,score_mismatch,score_match,score_mismatch,score_mismatch,0},
    // C
    {0,score_mismatch,score_mismatch,score_match,score_mismatch,0},
    // G
    {0,score_mismatch,score_mismatch,score_mismatch,score_match,0},
    // N
    {0,0,0,0,0,0}
};

__ocl_constant char endl = '\0';


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

__kernel void calc_alignment_smith_waterman(__global const char* read, __global const char* ref, __global char * alignment, __global short * alignment_index) {

	// Offset read and ref pointers
	read = read + g_id * read_length * VSIZE;
	ref = ref + g_id * ref_length * VSIZE;

	/*########################
	Init offsets and matrices
	##########################*/

	short16 best_read_pos = v_null;
	short16 best_ref_pos = v_null;

	short16 max_score = v_null;

	// Score matrix
	short16 score_matrix [(ref_length + 1) * SCORING_ROWS];

	// Backtracking matrix
	short16 backtrack_matrix [(ref_length + 1) * (read_length + 1)];

	// offset for first and second row
	int prev_row = 0;
	int cur_row = 1;
	int current_row_aln = 1;

	/*########################
	Calculate alignment and score matrices
	##########################*/

	for (short ref_pos = 0; ref_pos < ref_length + 1; ++ref_pos) {
		score_matrix[ref_pos] = v_null;
		backtrack_matrix[ref_pos] = v_start;
	}

	for (short read_pos = 0; read_pos < read_length; ++read_pos) {

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

		backtrack_matrix[current_row_aln * (ref_length + 1)] = v_start;
		score_matrix[cur_row * (ref_length + 1)] = v_null;

		//printf("Read base =%c and score=%i\n", read[read_pos],read_cache.s0);
		//printf("Read2 base =%c and score2=%i\n", read[read_pos+ read_length],read_cache.s1);

		for (short ref_pos = 0; ref_pos < ref_length; ++ref_pos) {

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

			short16 up = score_matrix[prev_row * (ref_length + 1) + ref_pos + 1] + v_score_gap_ref;
			short16 left = score_matrix[cur_row * (ref_length + 1) + ref_pos] + v_score_gap_read;
			short16 diag = score_matrix[prev_row * (ref_length + 1) + ref_pos] + base_cmp;

			short16 cur = max(up,max(left, max(diag, v_null)));

			//printf("%hi ",cur.s0);

			short16 pointer = v_start;

			// Reverse if then else if order for select
			pointer = select(pointer, v_left, cur == left);
			pointer = select(pointer, v_up, cur == up);
			pointer = select(pointer, v_diag, cur == diag);
			pointer = select(pointer, v_start, cur == v_null);

			backtrack_matrix[current_row_aln * (ref_length + 1) + ref_pos + 1] = pointer;

			//printf("%hi ", pointer.s7);

			best_read_pos = select(best_read_pos,(short16)(read_pos),cur > max_score);
			best_ref_pos = select(best_ref_pos,(short16)(ref_pos),cur > max_score);
			//printf("Ref pos: %hi, best ref pos: %hi, cur: %hi, max: %hi\n", ref_pos, best_ref_pos.s0, cur.s0, max_score.s0);

			max_score = max(max_score,cur);

			score_matrix[cur_row * (ref_length + 1) + ref_pos + 1] =  cur;

		}
		//printf("\n");

		//printf("\nMAX_SCORE = %hi\n",max_score.s0);
		prev_row = cur_row++;
		cur_row %= SCORING_ROWS;
		++current_row_aln;
	}

	/*########################
	Backtrack
	##########################*/

	short * v_read_pos = (short *)(&best_read_pos);
	short * v_ref_pos = (short *)(&best_ref_pos);

	//printf("Best read pos: %hi\n", best_read_pos.s7);
	//printf("Best ref pos: %hi\n", best_ref_pos.s7);

	for (int i = 0; i < VSIZE; ++i) {

		short read_pos = v_read_pos[i];
		short ref_pos = v_ref_pos[i];

		//printf("Best read pos %i : %hi\n", i, read_pos);
		//printf("Best ref pos %i : %hi\n", i, ref_pos);

		int aln_pos = aln_length - 2;

		short backtrack = ((short *)(&backtrack_matrix[(read_pos + 1) * (ref_length + 1) + ref_pos + 1]))[i];

		//printf("Backtrack %i : %hi\n", i, backtrack);

		while(backtrack != START) {

			//printf("Backtrack: %hi\n", backtrack);

			if (backtrack == UP) {
				alignment[(2 * aln_length * i) + aln_length + aln_pos] = '-';
				alignment[(2 * aln_length * i) + aln_pos] = read[read_pos + i * read_length];
				--read_pos;
			}
			if (backtrack == LEFT) {
				alignment[(2 * aln_length * i) + aln_pos] = '-';
				alignment[(2 * aln_length * i) + aln_length + aln_pos] = ref[ref_pos + i * ref_length];
				--ref_pos;
			}
			if (backtrack == DIAG) {
				alignment[(2 * aln_length * i) + aln_pos] = read[read_pos + i * read_length];
				alignment[ (2 * aln_length * i) + aln_length + aln_pos] = ref[ref_pos + i * ref_length];
				--read_pos;
				--ref_pos;
			}
			backtrack = ((short *)(&backtrack_matrix[(read_pos + 1) * (ref_length + 1) + ref_pos + 1]))[i];
			--aln_pos;
		}

		alignment_index[(2 * i)] = aln_pos + 1;
		alignment_index[(2 * i) + 1] = aln_pos + 1;
	}
}

__kernel void calc_alignment_needleman_wunsch(__global const char* read, __global const char* ref, __global char * alignment, __global short * alignment_index) {

	// Offset read and ref pointers
	read = read + g_id * read_length * VSIZE;
	ref = ref + g_id * ref_length * VSIZE;

	/*########################
	Init offsets and matrices
	##########################*/

	short16 best_read_pos = v_null;
	short16 best_ref_pos = v_null;

	short16 max_read_pos = (short16)(read_length - 1);
	short16 max_ref_pos = (short16)(ref_length - 1);

	short16 row_max = (short16)(SHORT_MIN);

	short16 global_row_max_index = (short16)(-1);

	short16 row_max_index = v_null;

	// Score matrix
	short16 score_matrix [(ref_length + 1) * SCORING_ROWS];

	// Backtracking matrix
	short16 backtrack_matrix [(ref_length + 1) * (read_length + 1)];

	// offset for first and second row
	int prev_row = 0;
	int cur_row = 1;
	int current_row_aln = 1;

	/*########################
	Calculate alignment and score matrices
	##########################*/

	for (short ref_pos = 0; ref_pos < ref_length + 1; ++ref_pos) {
		score_matrix[ref_pos] = v_null;
		backtrack_matrix[ref_pos] = v_start;
		//printf("%hi ",v_start.sF);
	}
	//printf("\n");

	for (short read_pos = 0; read_pos < read_length; ++read_pos) {

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

		backtrack_matrix[current_row_aln * (ref_length + 1)] = v_up;
		//printf("%hi ", backtrack_matrix[current_row_aln * (ref_length + 1)].sF);
		score_matrix[cur_row * (ref_length + 1)] = (short16)((read_pos + 1) * score_gap_ref);

		max_read_pos = select(max_read_pos, (short16)(read_pos - 1), (max_read_pos == (short16)(read_length - 1)) && convert_short16(read_cache == v_nullchar));
		global_row_max_index = select(global_row_max_index, row_max_index, max_read_pos + (short16)(1) == (short16)(read_pos));

		row_max = score_matrix[cur_row * (ref_length + 1)];
		row_max_index = v_null;

		//printf("Read base =%c and score=%i\n", read[read_pos],read_cache.s0);
		//printf("Read2 base =%c and score2=%i\n", read[read_pos+ read_length],read_cache.s1);

		for (short ref_pos = 0; ref_pos < ref_length; ++ref_pos) {

			char16 ref_cache = (char16)(char_to_score[ref[ref_pos]],
										char_to_score[ref[ref_pos + ref_length]],
										char_to_score[ref[ref_pos + 2 * ref_length]],
										char_to_score[ref[ref_pos + 3 * ref_length]],
										char_to_score[ref[ref_pos + 4 * ref_length]],
										char_to_score[ref[ref_pos + 5 * ref_length]],
										char_to_score[ref[ref_pos + 6 * ref_length]],
										char_to_score[ref[ref_pos + 7 * ref_length]],
										char_to_score[ref[ref_pos + 8 * ref_length]],
										char_to_score[ref[ref_pos + 9 * ref_length]],
										char_to_score[ref[ref_pos + 10 * ref_length]],
										char_to_score[ref[ref_pos + 11 * ref_length]],
										char_to_score[ref[ref_pos + 12 * ref_length]],
										char_to_score[ref[ref_pos + 13 * ref_length]],
										char_to_score[ref[ref_pos + 14 * ref_length]],
										char_to_score[ref[ref_pos + 15 * ref_length]]
			);

			short16 base_cmp = (short16)(base_score[read_cache.s0][ref_cache.s0],
					base_score[read_cache.s1][ref_cache.s1],
					base_score[read_cache.s2][ref_cache.s2],
					base_score[read_cache.s3][ref_cache.s3],
					base_score[read_cache.s4][ref_cache.s4],
					base_score[read_cache.s5][ref_cache.s5],
					base_score[read_cache.s6][ref_cache.s6],
					base_score[read_cache.s7][ref_cache.s7],
					base_score[read_cache.s8][ref_cache.s8],
					base_score[read_cache.s9][ref_cache.s9],
					base_score[read_cache.sA][ref_cache.sA],
					base_score[read_cache.sB][ref_cache.sB],
					base_score[read_cache.sC][ref_cache.sC],
					base_score[read_cache.sD][ref_cache.sD],
					base_score[read_cache.sE][ref_cache.sE],
					base_score[read_cache.sF][ref_cache.sF]
			);

			// up = matrix[prev_row * (ref_length + 1) + ref_pos + 1];
			// diag = matrix[prev_row * (ref_length + 1) + ref_pos];
			// left = matrix[current_row * (ref_length + 1) + ref_pos];

			short16 up = score_matrix[prev_row * (ref_length + 1) + ref_pos + 1] + v_score_gap_ref;
			short16 left = score_matrix[cur_row * (ref_length + 1) + ref_pos] + v_score_gap_read;
			short16 diag = score_matrix[prev_row * (ref_length + 1) + ref_pos] + base_cmp;

			short16 cur = max(up,max(left, diag));

			//printf("%hi ",cur.sF);

			short16 pointer = v_start;

			// Reverse if then else if order for select
			pointer = select(pointer, v_left, cur == left);
			pointer = select(pointer, v_up, cur == up);
			pointer = select(pointer, v_diag, cur == diag);

			backtrack_matrix[current_row_aln * (ref_length + 1) + ref_pos + 1] = pointer;

			max_ref_pos = select(max_ref_pos, (short16)(ref_pos - 1),(max_ref_pos == (ref_length - (short16)(1))) && convert_short16(ref_cache == v_nullchar));
			row_max_index = select(row_max_index, ref_pos, cur > row_max);
			row_max = select(row_max, cur, cur > row_max);

			score_matrix[cur_row * (ref_length + 1) + ref_pos + 1] =  cur;

		}
		//printf("\n");

		//printf("\nMAX_SCORE = %hi\n",max_score.s0);
		prev_row = cur_row++;
		cur_row %= SCORING_ROWS;
		++current_row_aln;
	}

	global_row_max_index = select(global_row_max_index, row_max_index, global_row_max_index < v_null);

	best_read_pos = max_read_pos;
	//best_ref_pos = select(global_row_max_index, max_ref_pos, max_ref_pos < global_row_max_index);
	best_ref_pos = select(global_row_max_index, max_ref_pos, max_ref_pos < global_row_max_index);

	/*########################
	Backtrack
	##########################*/

	short * v_read_pos = (short *)(&best_read_pos);
	short * v_ref_pos = (short *)(&best_ref_pos);

	//printf("Best read pos: %hi\n", best_read_pos.s7);
	//printf("Best ref pos: %hi\n", best_ref_pos.s7);

	for (int i = 0; i < VSIZE; ++i) {

		short read_pos = v_read_pos[i];
		short ref_pos = v_ref_pos[i];

		int aln_pos = aln_length - 2;

		short backtrack = ((short *)(&backtrack_matrix[(read_pos + 1) * (ref_length + 1) + ref_pos + 1]))[i];

		//printf("Backtrack %i : %hi\n", i, backtrack);

		while(backtrack != START) {

			//printf("Backtrack: %hi\n", backtrack);

			if (backtrack == UP) {
				alignment[(2 * aln_length * i) + aln_length + aln_pos] = '-';
				alignment[(2 * aln_length * i) + aln_pos] = read[read_pos + i * read_length];
				--read_pos;
			}
			if (backtrack == LEFT) {
				alignment[(2 * aln_length * i) + aln_pos] = '-';
				alignment[(2 * aln_length * i) + aln_length + aln_pos] = ref[ref_pos + i * ref_length];
				--ref_pos;
			}
			if (backtrack == DIAG) {
				alignment[(2 * aln_length * i) + aln_pos] = read[read_pos + i * read_length];
				alignment[ (2 * aln_length * i) + aln_length + aln_pos] = ref[ref_pos + i * ref_length];
				--read_pos;
				--ref_pos;
			}

			backtrack = ((short *)(&backtrack_matrix[(read_pos + 1) * (ref_length + 1) + ref_pos + 1]))[i];
			--aln_pos;
		}

		alignment_index[(2 * i)] = aln_pos + 1;
		alignment_index[(2 * i) + 1] = aln_pos + 1;
	}
}


__kernel void print_string(__global const char* A, __global const char* B){

	

	printf("Read %i is ",g_id);
	for (int i = 0; i < read_length; ++i) {
		printf("%c",A[g_id * read_length + i]);
	}
	
	printf("\nRef %i is ",g_id);
		for (int i = 0; i < ref_length; ++i) {
		printf("%c",B[g_id * ref_length + i]);
	}
}
