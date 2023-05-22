#ifndef KSW2_SSE_FUNCTIONS_H_
#define KSW2_SSE_FUNCTIONS_H_

void hw_set1_epi8(datatype *ris, datatype v)
{
#pragma HLS INLINE
set1_epi8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = v;
    }
}

void hw_blendv_epi8(datatype *ris, datatype *a, datatype *b, datatype *mask)
{
#pragma HLS INLINE
blendv_epi8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = mask[i] ? b[i] : a[i];
    }
}

void hw_mask_with_two(datatype *ris, datatype *a, datatype *mask)
{
#pragma HLS INLINE
mask_with_two:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = mask[i] ? 2 : a[i];
    }
}
void hw_cmpeq_epi8(datatype *ris, datatype *a, datatype *b)
{
#pragma HLS INLINE
cmpeq_epi8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = (a[i] == b[i]) ? 0xff : 0;
    }
}

void hw_cmpgt_epi8(datatype *ris, datatype *a, datatype *b)
{
#pragma HLS INLINE
cmpgt_epi8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = (a[i] > b[i]) ? 0xff : 0;
    }
}

void hw_cmpgt_epi8_zero(datatype *ris, datatype *a)
{
#pragma HLS INLINE
cmpgt_epi8_zero:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = (a[i] > 0) ? 0xff : 0;
    }
}

void hw_loadu_si128_8(datatype *ris, datatype *a)
{
#pragma HLS INLINE

    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i];
    }
}

void hw_loadu_si32_8(int32_t *ris, datatype *a)
{
#pragma HLS INLINE
loadu_si32_8:
    for (int i = 0; i < 4; i++)
    {
#pragma HLS UNROLL
        ris[i] = (int32_t)a[i];
    }
}

void hw_or_si128_8(datatype *ris, datatype *a, datatype *b)
{
#pragma HLS INLINE
or_si128_8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i] | b[i];
    }
}

void hw_and_si128_8(datatype *ris, datatype *a, datatype *b)
{
#pragma HLS INLINE
and_si128_8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i] & b[i];
    }
}

void hw_and_si128_8_1(datatype *ris, datatype *a)
{
#pragma HLS INLINE
and_si128_8_1:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i] & 1;
    }
}

void hw_and_si128_8_8(datatype *ris, datatype *a)
{
#pragma HLS INLINE
and_si128_8_8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i] & 8;
    }
}

void hw_and_si128_8_16(datatype *ris, datatype *a)
{
#pragma HLS INLINE
and_si128_8_16:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i] & 16;
    }
}

void hw_max_epi8(datatype *ris, datatype *a, datatype *b)
{
#pragma HLS INLINE
max_epi8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = (a[i] > b[i]) ? a[i] : b[i];
    }
}

void hw_min_epi8(datatype *ris, datatype *a, datatype *b)
{
#pragma HLS INLINE
min_epi8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = (a[i] < b[i]) ? a[i] : b[i];
    }
}

void hw_add_epi8(datatype *ris, datatype *a, datatype *b)
{
#pragma HLS INLINE
add_epi8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i] + b[i];
    }
}

void hw_sub_epi8(datatype *ris, datatype *a, datatype *b)
{
#pragma HLS INLINE
sub_epi8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i] - b[i];
    }
}

void hw_srli_si128_8(datatype *ris, datatype *a, datatype num)
{
#pragma HLS INLINE
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = 0;
    }
    for (int i = 0; i < 16 - num; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i + num];
    }
}

void hw_slli_si128_8(datatype *ris, datatype *a, datatype num)
{
#pragma HLS INLINE
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = 0;
    }
    for (int i = 0; i < 16 - num; i++)
    {
#pragma HLS UNROLL
        ris[i + num] = a[i];
    }
}

void hw_set1_epi32(int32_t *ris, int32_t a)
{ //
#pragma HLS INLINE
set1_epi32:
    for (int i = 0; i < 4; i++)
    {
#pragma HLS UNROLL
        ris[i] = a;
    }
}

void hw_setr_epi32(int32_t *ris, int32_t *a)
{
#pragma HLS INLINE
    ris[3] = a[0];
    ris[2] = a[1];
    ris[1] = a[2];
    ris[0] = a[3];
}

void hw_cvtsi32_si128(int32_t *ris, int32_t a)
{
#pragma HLS INLINE
    for (int i = 0; i < 4; i++)
    {
#pragma HLS UNROLL
        ris[i] = 0;
    }
    ris[0] = a;
}

void hw_loadu_si128_32(int32_t *ris, int32_t *a)
{ 
#pragma HLS INLINE
loadu_si128_32:
    for (int i = 0; i < 4; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i];
    }
}

void hw_add_epi32(int32_t *ris, int32_t *a, int32_t *b)
{
#pragma HLS INLINE
add_epi32:
    for (int i = 0; i < 4; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i] + b[i];
    }
}

void hw_add_epi32_8(int32_t *ris, int32_t *a, datatype *b)
{
#pragma HLS INLINE
add_epi32:
    for (int i = 0; i < 4; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i] + b[i];
    }
}

void hw_add_sub_epi32_32_8(int32_t *ris, int32_t *a, datatype *b)
{
#pragma HLS INLINE
add_epi32:
    for (int i = 0; i < 4; i++)
    {
#pragma HLS UNROLL
        ris[i] += -a[i] + b[i];
    }
}

void hw_add_epi32_set(int32_t *ris, int32_t *a, int32_t b)
{
#pragma HLS INLINE
add_epi32:
    for (int i = 0; i < 4; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i] + b;
    }
}

void hw_sub_epi32(int32_t *ris, int32_t *a, int32_t *b)
{
#pragma HLS INLINE
sub_epi32:
    for (int i = 0; i < 4; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i] - b[i];
    }
}

void hw_cmpgt_epi32(int32_t *ris, int32_t *a, int32_t *b)
{
#pragma HLS INLINE
cmpgt_epi32:
    for (int i = 0; i < 4; i++)
    {
#pragma HLS UNROLL
        if (a[i] > b[i])
            ris[i] = 0xffffffff; 
        else
            ris[i] = 0;
    }
}

void equal_hw_8(datatype *ris, datatype *a)
{
#pragma HLS INLINE
equal_8:
    for (int i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        ris[i] = a[i];
    }
}

void to_hw_8(datatype *ris, int32_t *a)
{
#pragma HLS INLINE
	ap_uint<32> tmp;
	tmp.range(7,0) = a[0];
	tmp.range(15,8) = a[1];
	tmp.range(23,16) = a[2];
	tmp.range(31,24) = a[3];
	*ris = (int32_t)tmp;
    //memcpy(ris, (datatype *)a, 4); // always copy 4 bytes
}

void hw_blendv_epi32(int32_t *ris, int32_t *a, int32_t *b, int32_t *mask)
{
#pragma HLS INLINE
blendv_epi_32:
    for (int i = 0; i < 4; i++)
    {
#pragma HLS UNROLL
        if (mask[i])
            ris[i] = b[i];
        else
            ris[i] = a[i];
    }
}

template <typename T>
//#ifdef __SYNTHESIS__
T to_matrix(T mat[MAX_SEQ_LEN>>4][1<<4], unsigned index)
//#else
//T to_matrix(T **mat, unsigned index)
//#endif
{
#pragma HLS INLINE
    return mat[(index) >> (4)][(index) & (15)];
}

int ksw_hw_apply_zdrop(ksw_hw_extz_t *ez, int is_rot, int32_t H, int a, int b, int zdrop, datatype e)
{
    int r, t;
    if (is_rot)
        r = a, t = b;
    else
        r = a + b, t = a;
    if (H > (int32_t)ez->max)
    {
        ez->max = H, ez->max_t = t, ez->max_q = r - t;
    }
    else if (t >= ez->max_t && r - t >= ez->max_q)
    {
        int tl = t - ez->max_t, ql = (r - t) - ez->max_q, l;
        l = tl > ql ? tl - ql : ql - tl;
        if (zdrop >= 0 && ez->max - H > zdrop + l * e)
        {
            ez->zdropped = 1;
            return 1;
        }
    }
    return 0;
}

#endif
