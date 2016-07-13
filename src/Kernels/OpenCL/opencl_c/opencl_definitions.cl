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
