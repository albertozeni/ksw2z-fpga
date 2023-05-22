#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "ap_int.h"
#include "datatypes.h"
#include "hls_stream.h"
#include "ksw2_sse_functions.h"

///////////////////////////////////////// da ksw2.h

#define HW_KSW_NEG_INF -0x40000000

#define HW_KSW_SCORE_ONLY 0x01	// don't record alignment path/cigar
#define HW_KSW_RIGHT 0x02		// right-align gaps
#define HW_KSW_GENERIC_SC 0x04	// without this flag: match/mismatch only; last symbol is a wildcard
#define HW_KSW_APPROX_MAX 0x08	// approximate max; this is faster with sse
#define HW_KSW_APPROX_DROP 0x10 // approximate Z-drop; faster with sse
#define HW_KSW_EXTZ_ONLY 0x40	// only perform extension
#define HW_KSW_REV_CIGAR 0x80	// reverse CIGAR in the output
#define HW_KSW_SPLICE_FOR 0x100
#define HW_KSW_SPLICE_REV 0x200
#define HW_KSW_SPLICE_FLANK 0x400

#define MAX_TLEN MAX_SEQ_LEN
#define MAX_QLEN MAX_SEQ_LEN

const unsigned int depth_stream = DEPTH_STREAM;
const unsigned int no_couples_per_stream = NO_COUPLES_PER_STREAM / NUM_CU;

void ksw_hw_reset_extz(ksw_hw_extz_t *ez)
{
#pragma HLS INLINE
	ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
	ez->max = 0, ez->score = ez->mqe = ez->mte = HW_KSW_NEG_INF;
	ez->m_cigar = 0, ez->n_cigar = 0, ez->zdropped = 0, ez->reach_end = 0;
}

void read(hls::stream<input_datatype> &input_ouput_stream, input_datatype *input_output_g)
{
#pragma HLS INLINE OFF
read_stream_loop:
	for (int i = 0; i < PACKAGES_PER_MAX_SEQ_LEN * 2; i++)
	{
#pragma HLS PIPELINE
		input_ouput_stream.write(input_output_g[i]);
	}
}

void read_wrapper(hls::stream<input_datatype> &input_ouput_stream, input_datatype *input_output_g, int n_pairs)
{
#pragma HLS PIPELINE off
read_loop:
	for (int i = 0; i < n_pairs + 1; i++, input_output_g += PACKAGES_PER_MAX_SEQ_LEN * 2) // one additional pair as we count the info of the scoring matrix
		read(input_ouput_stream, input_output_g);
}

void collector(hls::stream<int> score_local_stream[NUM_CU],
			   hls::stream<int> &final_score_stream, int num_couples)
{

loop_collector:
	for (int i = 0; i < num_couples / NUM_CU; i++)
	{
	loop_collector_inner:
		for (int j = 0; j < NUM_CU; j++)
		{
#pragma HLS PIPELINE
			int tmp = score_local_stream[j].read();
			final_score_stream.write(tmp);
		}
	}
}

void write(hls::stream<int> &final_score_stream, int n,
		   input_datatype *input_output, int n_pairs, int &to_send)
{

#pragma HLS INLINE OFF
	static int tmp[NUM_TMP_WRITE];
#pragma HLS BIND_STORAGE variable = tmp type = ram_1p impl = bram
	tmp[n & (NUM_TMP_WRITE - 1) /*%NUM_TMP_WRITE*/] = final_score_stream.read();
	if ((n > 0 && (((n + 1) & (NUM_TMP_WRITE - 1)) == 0)) || n == n_pairs - 1)
	{
		int iter = (to_send >= NUM_TMP_WRITE) ? NUM_TMP_WRITE : to_send;

		for (int i = 0; i < iter; i++)
		{
#pragma HLS pipeline
			//			printf("%d\n", n + i + 1 - iter);
			input_output[n + i + 1 - iter] = tmp[i];
		}

		to_send -= iter;
	}
}

void write_wrapper(hls::stream<int> &scores_stream, input_datatype *scores, int n_pairs)
{
#pragma HLS PIPELINE off
// #pragma HLS INLINE OFF
write_loop:
	int to_send = n_pairs;

loop_write_score_wrapper:
	for (int n = 0; n < n_pairs; n++)
	{
#pragma HLS PIPELINE off
		write(scores_stream, n, scores, n_pairs, to_send);
	}
	//	for (int i = 0; i < n_pairs; i++)
	//	{
	// #pragma HLS PIPELINE
	//		scores[i] = scores_stream.read();
	//	}
}

void dispatcher(hls::stream<input_datatype> &input_output_stream,
				hls::stream<input_datatype> local_input_output_stream[NUM_CU],
				int n_pairs)
{

	int idx = 0;

	input_datatype mat_local[PACKAGES_PER_MAX_SEQ_LEN * 2];

	for (int i = 0; i < PACKAGES_PER_MAX_SEQ_LEN * 2; i++)
		mat_local[i] = input_output_stream.read();

loop_matrix_broadcast:
	for (int i = 0; i < NUM_CU * PACKAGES_PER_MAX_SEQ_LEN * 2; i++)
	{
#pragma HLS PIPELINE
		local_input_output_stream[i / (PACKAGES_PER_MAX_SEQ_LEN * 2)].write(mat_local[i % (PACKAGES_PER_MAX_SEQ_LEN * 2)]);
	}

loop_dispatcher:
	for (int i = 0; i < n_pairs; i++, idx++)
	{
	loop_dispatcher_inner:
		for (int j = 0; j < PACKAGES_PER_MAX_SEQ_LEN * 2; j++)
		{
#pragma HLS PIPELINE
			input_datatype tmp_input = input_output_stream.read();
			local_input_output_stream[idx % NUM_CU].write(tmp_input);
		}
	}
}

void compute(int tlen, int qlen, hls::stream<input_datatype> &input_output_stream, int score, datatype *mat, datatype m,
			 datatype q, datatype e, int w, int zdrop, int end_bonus, int flag, ksw_hw_extz_t ez_local, int r, int t, int qe,
			 int n_col_, int tlen_, int qlen_, int last_st, int last_en, int wl, int wr, int max_sc, int min_sc,
			 int32_t H[MAX_TLEN >> 4][1 << 4], datatype *q_8, datatype *qe2_8, datatype *sc_mch_8, datatype *sc_mis_8, datatype *sc_N_8, datatype *m1_8, datatype *max_sc_8,
			 datatype u[MAX_TLEN >> 4][1 << 4], datatype v[MAX_TLEN >> 4][1 << 4], datatype x[MAX_TLEN][1 << 4], datatype y[MAX_TLEN][1 << 4], datatype s[MAX_TLEN >> 4][1 << 4], // ap_int<128> *p,
			 // int32_t *H132, int32_t *t_32, int32_t *tmp32, datatype *d8, datatype *d8p,
			 datatype *z8, datatype *z8p, datatype *xt18, datatype *xt18p, datatype *tmp8, datatype *tmp8p, datatype *vt18, datatype *vt18p, datatype *a8, datatype *ut8, datatype *b8, datatype *b8p,
			 int32_t *x1_32, int32_t *v1_32, datatype *x1_8, datatype *v1_8, datatype *sq8, datatype *st8, datatype *mask8, datatype *extra1, datatype *extra2, datatype zero,
			 int32_t *HH, int32_t *tt //, int32_t *max_H_32, int32_t *max_t_32, int32_t *qe_32//, int *off, int *off_end ,int32_t max_H, int32_t max_t, int32_t en1, int32_t i, int st, int en, int st0, int en0, int st_, int en_, datatype x1, datatype v1
			 ,
			 hls::stream<int> &scores_stream)
{
#pragma HLS INLINE OFF

#pragma HLS dependence variable = u inter false
#pragma HLS dependence variable = v inter false
#pragma HLS dependence variable = x inter false
#pragma HLS dependence variable = y inter false
#pragma HLS dependence variable = s inter false
// #pragma HLS dependence variable = p inter false
#pragma HLS dependence variable = H inter false

	//	read qr and sf from streams

	u_datatype qr[MAX_SEQ_LEN];
	u_datatype sf[MAX_SEQ_LEN];

#pragma HLS ARRAY_PARTITION variable = qr dim = 1 factor = 16 cyclic
#pragma HLS ARRAY_PARTITION variable = sf dim = 1 factor = 16 cyclic

#ifndef VITIS22
	input_datatype *qr_input = (input_datatype *)qr;
	input_datatype *sf_input = (input_datatype *)sf;
	for (int i = 0; i < PACKAGES_PER_MAX_SEQ_LEN; i++)
		qr_input[i] = input_output_stream.read();
	for (int i = 0; i < PACKAGES_PER_MAX_SEQ_LEN; i++)
		sf_input[i] = input_output_stream.read();
#else
	// ! for u55c (>= 2022.1)
	for (int i = 0; i < PACKAGES_PER_MAX_SEQ_LEN; i++)
	{
		input_datatype tmp_q = input_output_stream.read();
		input_datatype tmp_s = input_output_stream.read();
		for (int j = 0; j < CHAR_PER_PACKAGE; j++)
		{
#pragma HLS unroll
			qr[i * CHAR_PER_PACKAGE + j] = tmp_q.range(((j + 1) * CHAR_BITWIDTH) - 1, j * CHAR_BITWIDTH);
			sf[i * CHAR_PER_PACKAGE + j] = tmp_s.range(((j + 1) * CHAR_BITWIDTH) - 1, j * CHAR_BITWIDTH);
		}
	}
#endif
	ksw_hw_reset_extz(&ez_local);

	hw_set1_epi8(q_8, q);
	hw_set1_epi8(qe2_8, (q + e) * 2);
	hw_set1_epi8(sc_mch_8, mat[0]);
	hw_set1_epi8(sc_mis_8, mat[1]);
	mat[m * m - 1] == 0 ? hw_set1_epi8(sc_N_8, -e) : hw_set1_epi8(sc_N_8, mat[m * m - 1]);
	hw_set1_epi8(m1_8, m - 1);
	hw_set1_epi8(max_sc_8, (mat[0] + (q + e) * 2));

	if (w < 0)
		w = tlen > qlen ? tlen : qlen;
	wl = wr = w;
	tlen_ = (tlen + 15) / 16;
	n_col_ = qlen < tlen ? qlen : tlen;
	n_col_ = ((n_col_ < w + 1 ? n_col_ : w + 1) + 15) / 16 + 1;
	qlen_ = (qlen + 15) / 16;

matrix_scores_loop:
	for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t)
	{
		max_sc = max_sc > mat[t] ? max_sc : (int)mat[t];
		min_sc = min_sc < mat[t] ? min_sc : (int)mat[t];
	}

init_h:
	for (t = 0; t < MAX_TLEN / 16; t++)
	{
#pragma HLS PIPELINE
		for (int i = 0; i < 16; i++)
		{
#pragma HLS UNROLL
			H[t][i] = HW_KSW_NEG_INF;
		}
	}

// antidiags loop
antidiag_loop:
	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r)
	{
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		datatype x1, v1;
		u_datatype *qrr = qr + (qlen - 1 - r);

		// find the boundaries
		if (st < r - qlen + 1)
			st = r - qlen + 1;
		if (en > r)
			en = r;
		if (st < (r - wr + 1) >> 1)
			st = (r - wr + 1) >> 1; // take the ceil
		if (en > (r + wl) >> 1)
			en = (r + wl) >> 1; // take the floor
		if (st > en)
		{
			ez_local.zdropped = 1;
			break;
		}
		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
		// set boundary conditions
		if (st > 0)
		{
			if (st - 1 >= last_st && st - 1 <= last_en)
				x1 = to_matrix(x, st - 1), v1 = to_matrix(v, st - 1); // (r-1,s-1) calculated in the last round
			else
				x1 = v1 = 0; // not calculated; set to zeros
		}
		else
			x1 = 0, v1 = r ? q : zero;
		if (en >= r)
			y[r >> 4][r & 15] = 0, u[r >> 4][r & 15] = r ? q : zero;

		// loop fission: set scores first
		st_ = st / 16, en_ = en / 16;

		// core loop
		hw_cvtsi32_si128(x1_32, x1);
		hw_cvtsi32_si128(v1_32, v1);

		to_hw_8(x1_8, x1_32);
		to_hw_8(v1_8, v1_32);
		///////////////////////

		int32_t prev_H_en0 = to_matrix(H, en0);
		int32_t prev_H_en0_minus_1 = en0 >= 1 ? to_matrix(H, en0 - 1) : 0;

		for (int i = 0; i < 16; i++)
		{
#pragma HLS UNROLL
			HH[i] = KSW_NEG_INF;
		}

	core_loop:
		for (t = st_; t <= en_; t++)
		{
#pragma HLS PIPELINE
			for (int i = 0; i < 16; i++)
			{
#pragma HLS UNROLL
				tmp8[i] = mat[sf[(t << 4) + (st & 15) + i] * m + qrr[(t << 4) + (st & 15) + i]];
			}
			hw_loadu_si128_8(&s[t][st & 15], tmp8);

			//// __dp_block1

			hw_loadu_si128_8(z8p, s[t]); 
			hw_add_epi8(z8, z8p, qe2_8);
			hw_loadu_si128_8(xt18, x[t]);
			hw_srli_si128_8(tmp8, xt18, 15);
			hw_slli_si128_8(xt18p, xt18, 1);
			hw_or_si128_8(xt18, xt18p, x1_8);
			equal_hw_8(x1_8, tmp8);
			hw_loadu_si128_8(vt18, v[t]);
			hw_srli_si128_8(tmp8, vt18, 15);
			hw_slli_si128_8(vt18p, vt18, 1);
			hw_or_si128_8(vt18, vt18p, v1_8);
			equal_hw_8(v1_8, tmp8);
			hw_add_epi8(a8, xt18, vt18);
			hw_loadu_si128_8(ut8, u[t]);
			hw_loadu_si128_8(b8p, y[t]);
			hw_add_epi8(b8, b8p, ut8);

			/////////////////////

			hw_max_epi8(z8, z8, a8);

			/////////////////////////__dp_code_block2

			hw_max_epi8(z8, z8, b8);
			hw_min_epi8(z8, z8, max_sc_8);
			hw_sub_epi8(z8p, z8, vt18);
			hw_loadu_si128_8(u[t], z8p);
			hw_sub_epi8(z8p, z8, ut8);
			hw_loadu_si128_8(v[t], z8p);
			hw_sub_epi8(z8, z8, q_8);
			hw_sub_epi8(a8, a8, z8);
			hw_sub_epi8(b8, b8, z8);

			hw_cmpgt_epi8_zero(tmp8, a8);
			hw_and_si128_8(tmp8p, tmp8, a8);
			hw_loadu_si128_8(x[t], tmp8p);
			hw_cmpgt_epi8_zero(tmp8, b8);
			hw_and_si128_8(tmp8p, tmp8, b8);
			hw_loadu_si128_8(y[t], tmp8p);

			for (int i = 0; i < 16; i++)
			{
#pragma HLS UNROLL
				int idx = t * 16 + i;
				if (idx >= st0 && idx < en0)
				{
					int32_t tmp = H[t][i] + v[t][i] - qe - qe;
					H[t][i] = tmp;
					if (tmp > HH[i])
					{
						HH[i] = tmp;
						tt[i] = idx;
					}
				}
			}
		}

		int32_t max_H, max_t;

		if (r > 0)
		{
			int32_t u_en0 = to_matrix(u, en0);
			int32_t v_en0 = to_matrix(v, en0);
			max_t = en0;

			H[en0 / 16][en0 % 16] = max_H = en0 > 0 ? prev_H_en0_minus_1 + u_en0 : prev_H_en0 + v_en0; // special casing the last element

			for (int i = 0; i < 16; i++)
			{
#pragma HLS PIPELINE
				if (HH[i] > max_H)
				{
					max_H = HH[i];
					max_t = tt[i];
				}
			}
		}
		else
		{
			H[0][0] = v[0][0] - qe - qe;
			max_H = H[0][0];
			max_t = 0; // special casing r==0
		}

		int32_t H_en0 = to_matrix(H, en0);
		int32_t H_st0 = to_matrix(H, st0);

		if (en0 == tlen - 1 && H_en0 > ez_local.mte)
			ez_local.mte = H_en0, ez_local.mte_q = r - en;
		if (r - st0 == qlen - 1 && H_st0 > ez_local.mqe)
			ez_local.mqe = H_st0, ez_local.mqe_t = st0;
		if (ksw_hw_apply_zdrop(&ez_local, 1, max_H, r, max_t, zdrop, e))
		{
			break;
		}
		if (r == qlen + tlen - 2 && en0 == tlen - 1)
			score = to_matrix(H, tlen - 1);
		//			ez_local.score = H[tlen - 1];

		last_st = st, last_en = en;
		///////
	}

	scores_stream.write(score);
}

void compute_wrapper(int tlen, int qlen, hls::stream<input_datatype> &input_output_stream,
					 hls::stream<int> &scores_stream, int n_pairs, datatype q, datatype e, int w, int zdrop, int end_bonus, int flag, datatype m)
{

	int32_t H[MAX_TLEN >> 4][1 << 4];
	datatype u[MAX_TLEN >> 4][1 << 4];
	datatype v[MAX_TLEN >> 4][1 << 4];
	datatype x[MAX_TLEN >> 4][1 << 4];
	datatype y[MAX_TLEN >> 4][1 << 4];
	datatype s[MAX_TLEN >> 4][1 << 4];

	//	int32_t H132[4], t_32[4], tmp32[4];
	//	datatype d8[16], d8p[16];
	datatype z8[16], z8p[16], xt18[16], xt18p[16], tmp8[16], tmp8p[16], vt18[16], vt18p[16], a8[16], ut8[16], b8[16], b8p[16];
	ksw_hw_extz_t ez_local;
	//	int with_cigar = !(flag & HW_KSW_SCORE_ONLY), approx_max = !!(flag & HW_KSW_APPROX_MAX);
	datatype q_8[16], qe2_8[16], sc_mch_8[16], sc_mis_8[16], sc_N_8[16], m1_8[16], max_sc_8[16];
	int32_t x1_32[4], v1_32[4];
	datatype x1_8[16], v1_8[16];
	datatype sq8[16], st8[16], mask8[16];
	datatype extra1[16], extra2[16];
	datatype zero = 0;
	int32_t HH[16], tt[16];
	//	int32_t max_H_32[4], max_t_32[4], qe_32[4];
	//	int32_t max_H, max_t;
	int32_t en1, i;
	int st, en, st0, en0, st_, en_;
	datatype x1, v1;
	int score;
	int r, t, qe = q + e, n_col_, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc;
	//	int off[(MAX_QLEN + MAX_TLEN - 1) * 2];
	datatype mat[M * M];

#pragma HLS ARRAY_PARTITION variable = H complete dim = 2
// #pragma HLS BIND_STORAGE variable = H type = ram_t2p impl = uram
#pragma HLS BIND_STORAGE variable = u type = ram_t2p impl = uram
#pragma HLS BIND_STORAGE variable = v type = ram_t2p impl = uram
#pragma HLS BIND_STORAGE variable = x type = ram_t2p impl = uram
#pragma HLS BIND_STORAGE variable = y type = ram_t2p impl = uram
#pragma HLS BIND_STORAGE variable = s type = ram_t2p impl = uram
// #pragma HLS BIND_STORAGE variable = p type = ram_t2p impl = uram
// #pragma HLS BIND_STORAGE variable = off type = ram_t2p impl = uram

// #pragma HLS ARRAY_PARTITION variable = off dim = 1 factor = 16 cyclic
// #pragma HLS ARRAY_PARTITION variable=H dim=2 complete
// #pragma HLS ARRAY_RESHAPE variable = H dim = 2 complete
#pragma HLS ARRAY_RESHAPE variable = u dim = 2 complete
#pragma HLS ARRAY_RESHAPE variable = v dim = 2 complete
#pragma HLS ARRAY_RESHAPE variable = x dim = 2 complete
#pragma HLS ARRAY_RESHAPE variable = y dim = 2 complete
#pragma HLS ARRAY_RESHAPE variable = s dim = 2 complete

#pragma HLS ARRAY_PARTITION variable = z8 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = z8p dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = xt18 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = xt18p dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = tmp8 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = tmp8p dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = vt18 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = vt18p dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = a8 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = ut8 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = b8 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = b8p dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = x1_32 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = v1_32 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = x1_8 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = v1_8 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = sq8 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = st8 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = mask8 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = extra1 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = extra2 dim = 1 complete
// #pragma HLS ARRAY_PARTITION variable = max_H_32 dim = 1 complete
// #pragma HLS ARRAY_PARTITION variable = max_t_32 dim = 1 complete
// #pragma HLS ARRAY_PARTITION variable = qe_32 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = mat dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = HH dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = tt dim = 1 complete

#ifndef VITIS22
	// ! regular (<= 2020.2)
	input_datatype *mat_input = (input_datatype *)mat;
	// before computing read the scoring matrix
	for (int i = 0; i < PACKAGES_PER_MAX_SEQ_LEN * 2; i++)
	{
		if (i < 2)
			mat_input[i] = input_output_stream.read();
		else
			input_output_stream.read(); // empty garbage packets, might be used afterwards for additional info
	}
#else
	 ! for u55c (2022.1)
		for (int i = 0; i < PACKAGES_PER_MAX_SEQ_LEN; i++)
	 {
		 if (i < 2)
		 {
			 input_datatype tmp_m = input_output_stream.read();
			 for (int j = 0; j < CHAR_PER_PACKAGE; j++)
			 {
#pragma HLS unroll
				 mat[i * CHAR_PER_PACKAGE + j] = tmp_m.range(((j + 1) * CHAR_BITWIDTH) - 1, j * CHAR_BITWIDTH);
			 }
		 }
		 else
			 input_output_stream.read(); // empty garbage packets, might be used afterwards for additional info
		 input_output_stream.read();	 // empty garbage packets, might be used afterwards for additional info
	 }

//		 for (int i = 0; i < M * M; i++)
//		 {
//		 	printf("%d ", mat[i]);
//		 }
//		 printf("\n"); // DEBUG
#endif

compute_loop:
	for (int i = 0; i < n_pairs; i++)
	{
#pragma HLS PIPELINE off
		compute(tlen, qlen, input_output_stream, score, mat, m, q, e, w, zdrop, end_bonus, flag, ez_local, r, t,
				qe, n_col_, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc, H, q_8, qe2_8, sc_mch_8, sc_mis_8, sc_N_8, m1_8, max_sc_8, u, v, x, y, s, // p,
				// H132, t_32, tmp32, d8, d8p,
				z8, z8p, xt18, xt18p, tmp8, tmp8p, vt18, vt18p, a8, ut8, b8, b8p, x1_32, v1_32, x1_8, v1_8, sq8, st8, mask8,
				extra1, extra2, zero, HH, tt //, max_H_32, max_t_32, qe_32, off, off_end, max_H, max_t, en1, i, st, en, st0, en0, st_, en_, x1, v1
				,
				scores_stream);
	}
}

void ksw2_execute(hls::stream<input_datatype> &input_output_stream, input_datatype *input_output_g,
				  int n_pairs, int tlen, int qlen, hls::stream<input_datatype> local_input_output_stream[NUM_CU],
				  datatype q, datatype e, int w, int zdrop, int end_bonus, int flag, datatype m,
				  hls::stream<int> scores_local_stream[NUM_CU], hls::stream<int> &scores_stream)
{
#pragma HLS INLINE

	read_wrapper(input_output_stream, input_output_g, n_pairs);

	dispatcher(input_output_stream, local_input_output_stream, n_pairs);

	for (int i = 0; i < NUM_CU; i++)
	{
#pragma HLS unroll
		compute_wrapper(tlen, qlen, local_input_output_stream[i],
						//				, max_H, max_t, en1, i, st, en, st0, en0, st_, en_, x1, v1,
						scores_local_stream[i], n_pairs, q, e, w, zdrop, end_bonus, flag, m);
	}
	collector(scores_local_stream, scores_stream, n_pairs);

	write_wrapper(scores_stream, input_output_g, n_pairs);
}

// extern "C"
// {
void ksw2_extz2_hw(int qlen, int tlen, input_datatype *input_output_g, datatype m,				   //, const input_datatype *mat_global,
				   datatype q, datatype e, int w, int zdrop, int end_bonus, int flag, int n_pairs) // ksw_hw_extz_t *ez_out)
{
// depth is only for cosim purposes, query_l and target_l account for one additional couple for scoring infos
#pragma HLS INTERFACE m_axi offset = slave bundle = gmem0 port = input_output_g depth = 72

#pragma HLS INTERFACE s_axilite port = input_output_g bundle = control

#pragma HLS INTERFACE s_axilite port = qlen bundle = control
#pragma HLS INTERFACE s_axilite port = tlen bundle = control
#pragma HLS INTERFACE s_axilite port = m bundle = control
#pragma HLS INTERFACE s_axilite port = q bundle = control
#pragma HLS INTERFACE s_axilite port = e bundle = control
#pragma HLS INTERFACE s_axilite port = w bundle = control
#pragma HLS INTERFACE s_axilite port = zdrop bundle = control
#pragma HLS INTERFACE s_axilite port = end_bonus bundle = control
#pragma HLS INTERFACE s_axilite port = flag bundle = control
#pragma HLS INTERFACE s_axilite port = n_pairs bundle = control

#pragma HLS INTERFACE s_axilite port = return bundle = control
	printf("%d %d", depth_stream, no_couples_per_stream);
#pragma HLS dataflow

	static hls::stream<input_datatype> input_output_stream("input_stream");
#pragma HLS STREAM variable = input_output_stream depth = no_couples_per_stream
#pragma HLS BIND_STORAGE variable = input_output_stream type = fifo impl = bram

	static hls::stream<input_datatype> local_input_output_stream[NUM_CU];
#pragma HLS STREAM variable = local_input_output_stream depth = depth_stream
#pragma HLS BIND_STORAGE variable = local_input_output_stream type = fifo impl = bram

	static hls::stream<int> scores_stream("scores_stream");
#pragma HLS STREAM variable = scores_stream depth = depth_stream

	static hls::stream<int> scores_local_stream[NUM_CU];
#pragma HLS STREAM variable = scores_local_stream depth = no_couples_per_stream

	ksw2_execute(input_output_stream, input_output_g,
				 n_pairs, tlen, qlen, local_input_output_stream,
				 q, e, w, zdrop, end_bonus, flag, m, scores_local_stream, scores_stream);
}
// }
