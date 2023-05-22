#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>
#include <stdio.h>
#include "datatypes.h"
#include "xcl2.hpp"
#include "ksw2.h"

// #define NUM_KERNEL 10
// #define N_PAIRS 120

int NUM_KERNEL = 4;
int N_PAIRS = 0;

// HBM Banks requirements
#define MAX_DDR_BANKCOUNT 4
#define BANK_NAME(n) n | XCL_MEM_TOPOLOGY
const int bank[MAX_DDR_BANKCOUNT] = {
	BANK_NAME(0), BANK_NAME(1), BANK_NAME(2), BANK_NAME(3)};

#define NOW std::chrono::high_resolution_clock::now()

using namespace std;
using namespace chrono;

#define KSW_NEG_INF -0x40000000

#define KSW_EZ_SCORE_ONLY 0x01	// don't record alignment path/cigar
#define KSW_EZ_RIGHT 0x02		// right-align gaps
#define KSW_EZ_GENERIC_SC 0x04	// without this flag: match/mismatch only; last symbol is a wildcard
#define KSW_EZ_APPROX_MAX 0x08	// approximate max; this is faster with sse
#define KSW_EZ_APPROX_DROP 0x10 // approximate Z-drop; faster with sse
#define KSW_EZ_EXTZ_ONLY 0x40	// only perform extension
#define KSW_EZ_REV_CIGAR 0x80	// reverse CIGAR in the output
#define KSW_EZ_SPLICE_FOR 0x100
#define KSW_EZ_SPLICE_REV 0x200
#define KSW_EZ_SPLICE_FLANK 0x400

void ksw2_extz2_hw(int qlen, int tlen, input_datatype *input_output_g, datatype m, //, const input_datatype *mat_global,
				   datatype q, datatype e, int w, int zdrop, int end_bonus, int flag, int n_pairs);

void ksw2_execute_fpga(int qlen, int tlen, std::vector<uint8_t, aligned_allocator<uint8_t>> *input_output,
					   int8_t m, datatype q, datatype e, int w,
					   int zdrop, int end_bonus, int flag,
					   int n_pairs, std::string binary_file);

void copy_ez(ksw_hw_extz_t *ez_out, ksw_extz_t *ez)
{
	//#pragma HLS INLINE
	ez_out->max_q = ez->max_q;
	ez_out->max_t = ez->max_t;
	ez_out->mqe_t = ez->mqe_t;
	ez_out->mte_q = ez->mte_q;
	ez_out->max = ez->max;
	ez_out->score = ez->score;
	ez_out->mqe = ez->mqe;
	ez_out->mte = ez->mte;
	ez_out->m_cigar = ez->m_cigar;
	ez_out->n_cigar = ez->n_cigar;
	ez_out->zdropped = ez->zdropped;
	ez_out->reach_end = ez->reach_end;
}

//=======================================================================
//
// Common functions
//
//=======================================================================

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0 ? -a : a;
	b = b > 0 ? -b : b;
	for (i = 0; i < m - 1; ++i)
	{
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j ? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

static void global_aln_sw(const char *qseq_, const char *tseq_, int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
						  int w, int zdrop, int flag, ksw_extz_t *ez)
{
	int i, qlen, tlen;
	uint8_t *qseq, *tseq;
	ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
	ez->max = 0, ez->mqe = ez->mte = KSW_NEG_INF;
	ez->n_cigar = 0;
	qlen = MAX_SEQ_LEN; // strlen(qseq_);
	// printf("\ntarget: %s\n",qseq_); //
	tlen = MAX_SEQ_LEN; // strlen(tseq_);
	// printf("query: %s\n",tseq_); //
	qseq = (uint8_t *)calloc(qlen + 33, 1); // 32 for gaba
	tseq = (uint8_t *)calloc(tlen + 33, 1);
	for (i = 0; i < qlen; ++i)
		qseq[i] = seq_nt4_table[(uint8_t)qseq_[i]];
	for (i = 0; i < tlen; ++i)
		tseq[i] = seq_nt4_table[(uint8_t)tseq_[i]];

	/*//stampa delle sequenze
			printf("\n");
			for(i=0; i<qlen; i++)	printf("%d", qseq[i]);
			printf("\n");
			for(i=0; i<tlen; i++)	printf("%d", tseq[i]);
			printf("\n"); */
	flag &= KSW_EZ_SCORE_ONLY; // we are only computing the score
	// flag |= KSW_EZ_GENERIC_SC; //we are using the scoring matrix
	ksw_extz2_sse(0, qlen, (uint8_t *)qseq, tlen, (uint8_t *)tseq, m, mat, q, e, w, zdrop, 0, flag, ez);

	free(qseq);
	free(tseq);
}

//////hw
static void global_aln_hw(const char *qseq_, const char *tseq_, int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
						  int w, int zdrop, int flag, ksw_hw_extz_t *ez, std::vector<int, aligned_allocator<int>> *v_s, std::string binary_file)
{
	int i, qlen, tlen;

	qlen = (N_PAIRS)*MAX_SEQ_LEN;
	tlen = (N_PAIRS)*MAX_SEQ_LEN;

	std::vector<uint8_t, aligned_allocator<uint8_t>> totseq[NUM_KERNEL];
	for (int i = 0; i < NUM_KERNEL; i++)
		totseq[i].resize(qlen + tlen + MAX_SEQ_LEN * 2);
	for (unsigned n = 0; n < NUM_KERNEL; n++)
	{
		for (i = 0; i < m * m; i++)
		{
			totseq[n][i] = (datatype)mat[i];
		}
		for (i = 0; i < N_PAIRS; ++i)
		{
			for (int j = 0; j < MAX_SEQ_LEN; j++)
			{
				totseq[n][(1 + i) * MAX_SEQ_LEN * 2 + j] = seq_nt4_table[(uint8_t)qseq_[(i + 1) * MAX_SEQ_LEN - 1 - j]]; //(N_PAIRS*MAX_SEQ_LEN)-1-i]];
			}
			for (int j = 0; j < MAX_SEQ_LEN; j++)
			{
				totseq[n][(1 + i) * MAX_SEQ_LEN * 2 + MAX_SEQ_LEN + j] = seq_nt4_table[(uint8_t)tseq_[i * MAX_SEQ_LEN + j]];
			}
		}
	}
	flag |= KSW_EZ_SCORE_ONLY; // we are only computing the score
	flag |= KSW_EZ_GENERIC_SC; // we are using the scoring matrix
	ksw2_execute_fpga(MAX_SEQ_LEN, MAX_SEQ_LEN, totseq, m, q, e, w, zdrop, 0, flag, N_PAIRS, binary_file);
	for (int i = 0; i < NUM_KERNEL; i++)
	{
		uint8_t *totseq_array = totseq[i].data();
		for (int j = 0; j < N_PAIRS; j++)
			v_s[i][j] = (int)(((input_datatype *)totseq_array)[j]);
	}
	// free(qseq); free(tseq);
}

static void print_aln(const char *tname, const char *qname, ksw_extz_t *ez)
{
	printf("%s\t%s\t%d", tname, qname, ez->score);
	printf("\t%d\t%d\t%d", ez->max, ez->max_t, ez->max_q);
	if (ez->n_cigar > 0)
	{
		printf("\t%d", ez->n_cigar); //{
		int i;
		putchar('\t');
		for (i = 0; i < ez->n_cigar; ++i)
			printf("%d%c", ez->cigar[i] >> 4, "MID"[ez->cigar[i] & 0xf]);
	}
	putchar('\n');
}

/////////hw
static void print_aln_hw(const char *tname, const char *qname, ksw_hw_extz_t *ez)
{
	printf("%s\t%s\t%d", tname, qname, ez->score);
	printf("\t%d\t%d\t%d", ez->max, ez->max_t, ez->max_q);
	if (ez->n_cigar > 0)
	{
		printf("\t%d", ez->n_cigar);
		// int i;
		// putchar('\t');
		// 	for (i = 0; i < ez->n_cigar; ++i)
		// 		printf("%d%c", ez->cigar[i]>>4, "MID"[ez->cigar[i]&0xf]);
	}

	printf("\n");
}

static int check_res(int *comp_s, int *comp_t, int core_id)
{

	int flag = 0;
	for (int i = 0; !flag && i < N_PAIRS; i++)
	{

		if (comp_s[i] != comp_t[i])
		{
			printf("HW: %d, SW: %d, %d\n", comp_s[i], comp_t[i], i);
			flag = 1;
		}
	}

	if (flag == 1)
		printf("\nERROR ON CORE %d\n\n", core_id);

	return flag;
}

// static void random_generator(char* query, char* target, int qlen, int tlen)
// {
// 	char alphabet[4] = {'A', 'C', 'G', 'T'};
// 	int i=0;

// 	for(i=0; i<qlen; i++)	query[i] = alphabet[rand()%4];
// 	for(i=0; i<tlen; i++)	target[i] = alphabet[rand()%4];

// //	for(i=0; i<qlen; i++)	printf("%c",query[i]);
// //	putchar('\n');
// //	for(i=0; i<tlen; i++)	printf("%c", target[i]);
// //	putchar('\n');

// }

static void random_generator(char *query, char *target)
{
	char alphabet[5] = {'A', 'C', 'G', 'T', 'N'};
	for (int i = 0; i < (MAX_SEQ_LEN * N_PAIRS); i++)
		query[i] = alphabet[rand() % 4];
	for (int i = 0; i < (MAX_SEQ_LEN * N_PAIRS); i++)
		target[i] = alphabet[rand() % 4];
}

#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

#define EXPAND(a, m)                       \
	do                                     \
	{                                      \
		(m) = (m) ? (m) + ((m) >> 1) : 16; \
		REALLOC((a), (m));                 \
	} while (0)

typedef struct
{
	char *name;
	char *seq;
} named_seq_t;

//=======================================================================
//
// FPGA function calls
//
//=======================================================================

void ksw2_execute_fpga(int qlen, int tlen, std::vector<uint8_t, aligned_allocator<uint8_t>> *input_output,
					   int8_t m, datatype q, datatype e, int w,
					   int zdrop, int end_bonus, int flag,
					   int n_pairs, std::string binary_file)
{

	cl_int err;
	cl::CommandQueue queue;
	std::string krnl_name = "ksw2_extz2_hw";
	std::vector<cl::Kernel> krnls(NUM_KERNEL);
	cl::Context context;
	auto devices = xcl::get_xil_devices();

	// read_binary_file() command will find the OpenCL binary file created using the
	// V++ compiler load into OpenCL Binary and return pointer to file buffer.
	auto fileBuf = xcl::read_binary_file(binary_file);

	cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
	int valid_device = 0;
	for (unsigned int i = 0; i < devices.size(); i++)
	{
		auto device = devices[i];
		// Creating Context and Command Queue for selected Device
		OCL_CHECK(err, context = cl::Context(device, NULL, NULL, NULL, &err));
		OCL_CHECK(err,
				  queue = cl::CommandQueue(context,
										   device,
										   CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE |
											   CL_QUEUE_PROFILING_ENABLE,
										   &err));

		std::cout << "Trying to program device[" << i
				  << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
		cl::Program program(context, {device}, bins, NULL, &err);
		if (err != CL_SUCCESS)
		{
			std::cout << "Failed to program device[" << i
					  << "] with xclbin file!\n";
		}
		else
		{
			std::cout << "Device[" << i << "]: program successful!\n";
			//  Creating Kernel object using Compute unit names

			for (int i = 0; i < NUM_KERNEL; i++)
			{
				std::string cu_id = std::to_string(i + 1);
				std::string krnl_name_full =
					krnl_name + ":{" + "ksw2_extz2_hw_" + cu_id + "}";

				printf("Creating a kernel [%s] for CU(%d)\n",
					   krnl_name_full.c_str(),
					   i + 1);

				// Here Kernel object is created by specifying kernel name along with compute unit.
				// For such case, this kernel object can only access the specific Compute unit

				OCL_CHECK(err,
						  krnls[i] = cl::Kernel(
							  program, krnl_name_full.c_str(), &err));
			}
			valid_device++;
			break; // we break because we found a valid device
		}
	}
	if (valid_device == 0)
	{
		std::cout << "Failed to program any device found, exit!\n";
		exit(EXIT_FAILURE);
	}

	// SETUP DATA TRANSFER
	std::vector<cl_mem_ext_ptr_t> buffer_input_output_ext(NUM_KERNEL);
	std::vector<cl::Buffer> buffer_input_output(NUM_KERNEL);

	for (int i = 0; i < NUM_KERNEL; i++)
	{
		buffer_input_output_ext[i].obj = input_output[i].data();
		buffer_input_output_ext[i].param = 0;
		buffer_input_output_ext[i].flags = bank[i];
	}

	// These commands will allocate memory on the FPGA. The cl::Buffer objects can
	// be used to reference the memory locations on the device.
	// Creating Buffers
	for (int i = 0; i < NUM_KERNEL; i++)
	{
		OCL_CHECK(err,
				  buffer_input_output[i] =
					  cl::Buffer(context,
								 CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
									 CL_MEM_USE_HOST_PTR,
								 (N_PAIRS+1)*2* MAX_SEQ_LEN,
								 &buffer_input_output_ext[i],
								 &err));
	}

	// Copy input data to Device Global Memory
	for (int i = 0; i < NUM_KERNEL; i++)
	{
		OCL_CHECK(err,
				  err = queue.enqueueMigrateMemObjects(
					  {buffer_input_output[i]},
					  0 /* 0 means from host*/));
	}
	queue.finish();
	double kernel_time_in_sec = 0, result = 0;
	std::chrono::duration<double> kernel_time(0);
	auto kernel_start = NOW;

	for (int i = 0; i < NUM_KERNEL; i++)
	{
		// Setting kernel arguments

		OCL_CHECK(err, err = krnls[i].setArg(0, qlen));
		OCL_CHECK(err, err = krnls[i].setArg(1, tlen));
		OCL_CHECK(err, err = krnls[i].setArg(2, buffer_input_output[i]));
		OCL_CHECK(err, err = krnls[i].setArg(3, m));
		OCL_CHECK(err, err = krnls[i].setArg(4, q));
		OCL_CHECK(err, err = krnls[i].setArg(5, e));
		OCL_CHECK(err, err = krnls[i].setArg(6, w));
		OCL_CHECK(err, err = krnls[i].setArg(7, zdrop));
		OCL_CHECK(err, err = krnls[i].setArg(8, end_bonus));
		OCL_CHECK(err, err = krnls[i].setArg(9, flag));
		OCL_CHECK(err, err = krnls[i].setArg(10, n_pairs));

		// Invoking the kernel
		OCL_CHECK(err, err = queue.enqueueTask(krnls[i]));
	}

	queue.finish();

	double TOT_CELLS_HW = (MAX_SEQ_LEN * MAX_SEQ_LEN);
	TOT_CELLS_HW *= N_PAIRS;
	TOT_CELLS_HW *= NUM_KERNEL;
	auto kernel_end = NOW;
	kernel_time = std::chrono::duration<double>(kernel_end - kernel_start);
	kernel_time_in_sec = kernel_time.count();
	std::cerr << "\nFPGA TIME: " << kernel_time_in_sec << "s" << std::endl;
	std::cerr << "GCUPS FPGA: " << ((TOT_CELLS_HW / kernel_time_in_sec) / 1E9) << std::endl;

	for (int i = 0; i < NUM_KERNEL; i++)
	{
		OCL_CHECK(err,
				  err = queue.enqueueMigrateMemObjects(
					  {buffer_input_output[i]},
					  CL_MIGRATE_MEM_OBJECT_HOST));
	}
	queue.finish();

	/////////////////////////////////////////////////////////////////////////
	// printf("%d\n", v_s[0]);
	// //host mem reset
	// free(scoreLeft);
	// free(scoreRight);
}

//=======================================================================
//
// Function call main
//
//=======================================================================

int main(int argc, char *argv[])
{
	if(argc != 5){
		fprintf(stderr, "Usage: ./host-ksw2 bitstream npairs w zdrop\n");
		exit(-1);
	}
	string binary_file = argv[1];
	std::cerr << "[KSW2::] Read binary file: " << binary_file << std::endl;
	int8_t a = 2, b = 4, q = 4, e = 2, q2 = 13, e2 = 1;
	int c, i, pair = 1, w = -1, flag = 0, rep = 1, zdrop = -1, no_kalloc = 0;
	N_PAIRS = atoi(argv[2]);
	w = atoi(argv[3]);
	zdrop = atoi(argv[4]);
	char *s;
	int8_t mat[25];
	// ksw_extz_t ez;
	ksw_hw_extz_t ez_hw;

	// comparazione
	int comp_t[4];
	int comp_s[4];

	// generazione random

	int seed = time(NULL);
	printf("Seed: %d\n", seed);
	srand(seed);

	int random_qlen, random_tlen;
	random_qlen = rand() % MAX_SEQ_LEN;
	random_tlen = rand() % MAX_SEQ_LEN;

	// char random_query[MAX_SEQ_LEN * N_PAIRS];
	// char random_target[MAX_SEQ_LEN * N_PAIRS];

	// char query_test[MAX_SEQ_LEN];
	// char target_test[MAX_SEQ_LEN];

	char *random_query, *random_target;

	random_query = (char *)malloc(sizeof(char) * MAX_SEQ_LEN * N_PAIRS);
	random_target = (char *)malloc(sizeof(char) * MAX_SEQ_LEN * N_PAIRS);

	random_generator(random_query, random_target);

	memset(&ez_hw, 0, sizeof(ksw_hw_extz_t));
	ksw_gen_simple_mat(5, mat, a, -b);

	char *qpoint, *tpoint;

	qpoint = &random_query[0];
	tpoint = &random_target[0];

	// ksw_extz_t ris_test[N_PAIRS];
	// ksw_hw_extz_t ris_src[N_PAIRS];

	// int v_s[N_PAIRS];
	std::vector<int, aligned_allocator<int>> v_s[NUM_KERNEL];
	for (int i = 0; i < NUM_KERNEL; i++)
		v_s[i].resize(N_PAIRS);
	std::vector<int, aligned_allocator<int>> v_t(N_PAIRS);
	// int v_t[N_PAIRS];

	// copy_ez(&ez_hw, &ez);
	// ez_hw.score = 100;

	global_aln_hw(random_query, random_target, 5, mat, q, e, q2, e2, w, zdrop, flag, &ez_hw, v_s, binary_file);

	auto sw_start = NOW;

	// #pragma omp parallel for
	for (int j = 0; j < N_PAIRS; j++)
	{
		char *query_test, *target_test;

		query_test = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
		target_test = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
		ksw_extz_t ez;
		memset(&ez, 0, sizeof(ksw_extz_t));
		memcpy(query_test, qpoint + MAX_SEQ_LEN * j, sizeof(char) * MAX_SEQ_LEN);
		memcpy(target_test, tpoint + MAX_SEQ_LEN * j, sizeof(char) * MAX_SEQ_LEN);

		global_aln_sw(query_test, target_test, 5, mat, q, e, q2, e2, w, zdrop, flag, &ez);

		// ris_test[j].max_q = ez.max_q;
		// ris_test[j].max_t = ez.max_t;
		// ris_test[j].mqe_t = ez.mqe_t;
		// ris_test[j].mte_q = ez.mte_q;
		// ris_test[j].max = ez.max;
		// ris_test[j].score = ez.score;
		// ris_test[j].mqe = ez.mqe;
		// ris_test[j].mte = ez.mte;
		// ris_test[j].m_cigar = ez.m_cigar;
		// ris_test[j].n_cigar = ez.n_cigar;
		// ris_test[j].zdropped = ez.zdropped;
		// ris_test[j].reach_end = ez.reach_end;

		// print_aln("test", "", &ez);

		v_t[j] = ez.score;
		free(query_test);
		free(target_test);
	}

	double TOT_CELLS_SW = (MAX_SEQ_LEN * MAX_SEQ_LEN);
	TOT_CELLS_SW *= N_PAIRS;
	auto sw_end = NOW;
	std::chrono::duration<double> sw_time = std::chrono::duration<double>(sw_end - sw_start);
	double sw_time_in_sec = sw_time.count();
	std::cerr << "CPU TIME: " << sw_time_in_sec << "s" << std::endl;
	std::cerr << "GCUPS CPU: " << ((TOT_CELLS_SW / sw_time_in_sec) / 1E9) << std::endl;

	int error = 0;
	for (int i = 0; !error && i < NUM_KERNEL; i++)
		error |= check_res(v_s[i].data(), v_t.data(), i);
	if (!error)
		std::cerr << "ALL RESULTS OK" << std::endl;
	return 0;
}
