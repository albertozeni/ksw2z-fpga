#ifndef DATATYPES_H
#define DATATYPES_H

#include "ap_int.h"

// ! tune these for performance
#define NUM_TMP_WRITE 256
#define NUM_CU 14
#define MAX_SEQ_LEN 256
#define DEPTH_STREAM 100

#define LOG_2_PORT_WIDTH 8
#define LOG_2_BITS_PER_CHAR 3
#define CHAR_BITWIDTH 8
#define CHAR_PER_PACKAGE ((1 << LOG_2_PORT_WIDTH) >> 3)
#define PACKAGES_PER_MAX_SEQ_LEN (MAX_SEQ_LEN / CHAR_PER_PACKAGE)
#define NO_COUPLES_PER_STREAM ((DEPTH_STREAM / (PACKAGES_PER_MAX_SEQ_LEN * 2)) * NUM_CU)


#define KSW_NEG_INF -0x40000000
#define M 8 // dimension of the scoring matrix is M*M

typedef uint8_t u_datatype;
typedef int8_t datatype;
typedef ap_uint<1 << LOG_2_PORT_WIDTH> input_datatype;

typedef struct
{ //
	uint32_t max;
	uint32_t zdropped;
	int max_q, max_t; // max extension coordinate
	int mqe, mqe_t;	  // max score when reaching the end of query
	int mte, mte_q;	  // max score when reaching the end of target
	int score;		  // max score reaching both ends; may be HW_KSW_NEG_INF
	int m_cigar, n_cigar;
	int reach_end;
	int pat[4];
	//	uint32_t cigar[MAX_TLEN+MAX_QLEN];
} ksw_hw_extz_t;

#endif
