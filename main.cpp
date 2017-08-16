#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "Fourier.h"
#include "complex.h"
#include "line_fitting.h"
#include "viterbi_decoder.h"
#include "crc32_80211.h"
#include "common.h"
#include "mapper.h"
#include "sse_mathfun.h"

#ifdef WIN32
#include <Windows.h>
#include <intrin.h>
#include <tmmintrin.h>
#else
#include <x86intrin.h>
#include <tmmintrin.h>
#endif

using namespace gr::ieee802_11;

#define countof(x) (sizeof(x)/sizeof(x[0]))

const int CORDIC_TBL_COUNT = 100000;
const int MAX_SYMBOLS_PER_SAMPLE = 200;
complex cordic_tbl[CORDIC_TBL_COUNT];
int16_t cordic_tbl16[CORDIC_TBL_COUNT*2];
int interleave_pattern[4][288];	// [BPSK, QPSK, 16QAM, 64QAM] [max: 6*48bits]

int init_cordic()
{
	for(int i=0; i<CORDIC_TBL_COUNT; i++)
	{
		float phase = i*2*PI/CORDIC_TBL_COUNT;
		cordic_tbl[i] = complex(cos(phase), sin(phase));
		cordic_tbl16[i*2] = cordic_tbl[i].real * 32767;
		cordic_tbl16[i*2+1] = cordic_tbl[i].image * 32767;
	}
	return 0;
}

complex & cordic(float phase)
{
	static const float pi2 = acos(-1.0)*2;
	while(phase<0)
		phase += pi2;
	while(phase>pi2)
		phase -= pi2;
	int idx = int(phase/pi2 * (CORDIC_TBL_COUNT-1));
	assert(idx>=0 && idx <CORDIC_TBL_COUNT);
	return cordic_tbl[idx];
}

int16_t* cordic16(float phase)
{
	static const float pi2 = acos(-1.0)*2;
	while(phase<0)
		phase += pi2;
	while(phase>pi2)
		phase -= pi2;
	int idx = int(phase/pi2 * (CORDIC_TBL_COUNT-1));
	assert(idx>=0 && idx <CORDIC_TBL_COUNT);
	return cordic_tbl16+idx*2;
}


float correalation(float *s1, float *s2, int count)
{
	__m128 sum_a0 = _mm_setzero_ps();
	__m128 sum_a1 = _mm_setzero_ps();
	__m128 sum_p = _mm_setzero_ps();
	for(int i=0; i<count; i+=8)
	{
		__m128 a1 = _mm_loadu_ps(s1+i);
		__m128 a2 = _mm_loadu_ps(s1+4+i);
		__m128 b1 = _mm_loadu_ps(s2+i);
		__m128 b2 = _mm_loadu_ps(s2+4+i);

		__m128 tmp1 = _mm_mul_ps(a1, b1);
		__m128 tmp2 = _mm_mul_ps(a2, b2);
		tmp1 = _mm_hadd_ps(tmp1, tmp2);
		sum_a0 = _mm_add_ps(sum_a0, tmp1);

		tmp1 = _mm_mul_ps(a1, a1);
		sum_p = _mm_add_ps(sum_p, tmp1);
		tmp1 = _mm_mul_ps(a2, a2);
		sum_p = _mm_add_ps(sum_p, tmp1);

		tmp1 = _mm_shuffle_ps(b1, b1, _MM_SHUFFLE(2,3,0,1));
		tmp2 = _mm_shuffle_ps(b2, b2, _MM_SHUFFLE(2,3,0,1));

		tmp1 = _mm_mul_ps(a1, tmp1);
		tmp2 = _mm_mul_ps(a2, tmp2);
		tmp1 = _mm_hsub_ps(tmp1, tmp2);
		sum_a1 = _mm_add_ps(sum_a1, tmp1);
	}

	sum_a0 = _mm_hadd_ps(sum_a0, sum_a1);
	sum_a0 = _mm_hadd_ps(sum_a0, sum_a0);
	sum_a0 = _mm_mul_ps(sum_a0, sum_a0);
	sum_a0 = _mm_hadd_ps(sum_a0, sum_a0);
	sum_a0 = _mm_sqrt_ps(sum_a0);
	sum_p = _mm_hadd_ps(sum_p, sum_p);
	sum_p = _mm_hadd_ps(sum_p, sum_p);
	sum_a0 = _mm_div_ps(sum_a0, sum_p);


	return *(float*)&sum_a0;
}


// abj = a * b.conjugate()
// amp = |a|^2, no image part.
// count: number of int16_t, not number of complex number
int complex16_abj_amp(int16_t *a, int16_t *b, int16_t *abj, int16_t *amp, int count)
{
	__m128i mff = _mm_set1_epi16(0xff);
	__m128i m_min = _mm_set1_epi16(0);

	for(int i=0; i<count; i+=16)
	{
		__m128i * pa = (__m128i*)(a+i);
		__m128i * pb = (__m128i*)(b+i);
		__m128i * pout = (__m128i*)(abj+i);

		__m128i ma1 = _mm_loadu_si128(pa);
		__m128i ma2 = _mm_loadu_si128(pa+1);
		__m128i mb1 = _mm_loadu_si128(pb);
		__m128i mb2 = _mm_loadu_si128(pb+1);
		ma1 = _mm_srai_epi16(ma1, 8);
		ma2 = _mm_srai_epi16(ma2, 8);
		mb1 = _mm_srai_epi16(mb1, 8);
		mb2 = _mm_srai_epi16(mb2, 8);

		__m128i tmp1 = _mm_mullo_epi16(ma1, ma1);
		__m128i tmp2 = _mm_mullo_epi16(ma2, ma2);		// RaRa IaIa ...
		tmp2 = _mm_hadd_epi16(tmp1, tmp2);				// RaRa+IaIa ...
		tmp2 = _mm_max_epi16(m_min, tmp2);				// at least 20*20
		_mm_storeu_si128((__m128i*)(amp+i/2), tmp2);

		tmp1 = _mm_mullo_epi16(ma1, mb1);
		tmp2 = _mm_mullo_epi16(ma2, mb2);				// RaRb IaIb ...

		__m128i mreal = _mm_hadd_epi16(tmp1, tmp2);		// RaRb+IaIb ...

		tmp1 = _mm_srli_epi32(ma1, 16);
		ma1 = _mm_slli_epi32(ma1, 16);
		ma1 = _mm_xor_si128(ma1, tmp1);

		tmp1 = _mm_srli_epi32(ma2, 16);
		ma2 = _mm_slli_epi32(ma2, 16);
		ma2 = _mm_xor_si128(ma2, tmp1);					// swap real and image of a, -->> Ia Ra ...

		tmp1 = _mm_mullo_epi16(ma1, mb1);
		tmp2 = _mm_mullo_epi16(ma2, mb2);				// IaRb RaIb ...
		__m128i mimage = _mm_hsub_epi16(tmp1, tmp2);	// IaRb-RaIb ...

		tmp1 = _mm_unpacklo_epi16(mreal, mimage);
		tmp2 = _mm_unpackhi_epi16(mreal, mimage);		// interleave

		_mm_storeu_si128(pout, tmp1);
		_mm_storeu_si128(pout+1, tmp2);
	}

	return count;
}

int complex16_frequency_offset(int16_t *samples, int16_t *out, int count, float offset)
{
	for(int i=0; i<count; i+=16)
	{
		__m128 f1 = _mm_set_ps(i/2+0, i/2+1, i/2+2, i/2+3);
		__m128 f2 = _mm_set1_ps(4);
		__m128 fre = _mm_set1_ps(offset);
		f2 = _mm_add_ps(f1, f2);
		f1 = _mm_mul_ps(f1, fre);
		f2 = _mm_mul_ps(f2, fre);
		__m128 s1,s2,c1,c2;
		sincos_ps(f1, &s1, &c1);
		sincos_ps(f2, &s2, &c2);
		f1 = _mm_set1_ps(32767);
		s1 = _mm_mul_ps(f1, s1);
		s2 = _mm_mul_ps(f1, s2);
		c1 = _mm_mul_ps(f1, c1);
		c2= _mm_mul_ps(f1, c2);

		__m128i s1i = _mm_cvtps_epi32(s1);
		__m128i s2i = _mm_cvtps_epi32(s2);
		__m128i c1i = _mm_cvtps_epi32(c1);
		__m128i c2i = _mm_cvtps_epi32(c2);

		s1i = _mm_packs_epi32(s1i, s2i);
		c1i = _mm_packs_epi32(c1i, c2i);


		__m128i mb1 = _mm_unpacklo_epi16(c1i, s1i);
		__m128i mb2 = _mm_unpackhi_epi16(c1i, s1i);


		__m128i ma1 = _mm_loadu_si128((__m128i*)(samples+i));
		__m128i ma2 = _mm_loadu_si128((__m128i*)(samples+i+8));

		__m128i tmp1 = _mm_mulhi_epi16(ma1, mb1);
		__m128i tmp2 = _mm_mulhi_epi16(ma2, mb2);			// RaRb IaIb ...

		__m128i mreal = _mm_hsub_epi16(tmp1, tmp2);			// RaRb-IaIb ...

		tmp1 = _mm_srli_epi32(ma1, 16);
		ma1 = _mm_slli_epi32(ma1, 16);
		ma1 = _mm_xor_si128(ma1, tmp1);

		tmp1 = _mm_srli_epi32(ma2, 16);
		ma2 = _mm_slli_epi32(ma2, 16);
		ma2 = _mm_xor_si128(ma2, tmp1);					// swap real and image of a, -->> Ia Ra ...

		tmp1 = _mm_mulhi_epi16(ma1, mb1);
		tmp2 = _mm_mulhi_epi16(ma2, mb2);				// IaRb RaIb ...
		__m128i mimage = _mm_hadd_epi16(tmp1, tmp2);	// IaRb+RaIb ...

		tmp1 = _mm_unpacklo_epi16(mreal, mimage);
		tmp2 = _mm_unpackhi_epi16(mreal, mimage);		// interleave

		_mm_storeu_si128((__m128i*)(out+i), tmp1);
		_mm_storeu_si128((__m128i*)(out+i+8), tmp2);
	}

	return 0;
}


int complex_frequency_offset(complex *samples, int count, float offset, float phase = 0, complex *out = NULL)
{
	if (out == NULL)
		out = samples;

	__m128 mphase = _mm_set1_ps(phase);
	__m128 fre = _mm_set1_ps(offset);

	for(int i=0; i<count; i+=4)
	{
		__m128 f1 = _mm_set_ps(i+3, i+2, i+1, i+0);
		f1 = _mm_mul_ps(f1, fre);
		f1 = _mm_add_ps(f1, mphase);
		__m128 s1,c1;
		sincos_ps(f1, &s1, &c1);

		__m128 mb1 = _mm_unpacklo_ps(c1, s1);
		__m128 mb2 = _mm_unpackhi_ps(c1, s1);

		__m128 ma1 = _mm_loadu_ps((float*)(samples+i));
		__m128 ma2 = _mm_loadu_ps((float*)(samples+i+2));

		__m128 tmp1 = _mm_mul_ps(ma1, mb1);
		__m128 tmp2 = _mm_mul_ps(ma2, mb2);			// RaRb IaIb ...

		__m128 mreal = _mm_hsub_ps(tmp1, tmp2);		// RaRb-IaIb ...


		mb1 = _mm_unpacklo_ps(s1, c1);
		mb2 = _mm_unpackhi_ps(s1, c1);

		tmp1 = _mm_mul_ps(ma1, mb1);
		tmp2 = _mm_mul_ps(ma2, mb2);				// RaIb IaRb ...
		__m128 mimage = _mm_hadd_ps(tmp1, tmp2);	// RaIb+IaRb ...

		tmp1 = _mm_unpacklo_ps(mreal, mimage);
		tmp2 = _mm_unpackhi_ps(mreal, mimage);		// interleave

		_mm_storeu_ps((float*)(out+i), tmp1);
		_mm_storeu_ps((float*)(out+i+2), tmp2);
	}

	return 0;
}

// count: number of complex
int complex_from_int16(complex *out, int16_t *in, int count, bool swap = false)
{
	__m128i zero = _mm_setzero_si128();
	for(int i=0; i<count*2; i+=8)
	{
		__m128i m = _mm_loadu_si128((__m128i*)(in+i));
		if (swap)
		{
			__m128i a = _mm_srli_epi32(m, 16);
			__m128i b = _mm_slli_epi32(m, 16);
			m = _mm_or_si128(a, b);
		}
		__m128i sign = _mm_cmplt_epi16(m, zero);
		__m128i m1 = _mm_unpacklo_epi16(m, sign);
		__m128i m2 = _mm_unpackhi_epi16(m, sign);

		__m128 m1f = _mm_cvtepi32_ps(m1);
		__m128 m2f = _mm_cvtepi32_ps(m2);

		_mm_storeu_ps((float*)(out+i/2), m1f);
		_mm_storeu_ps((float*)(out+i/2+2), m2f);
	}

	return 0;
}

int init_interleaver_pattern()
{
	int bpsc_tbl[4] = {1,2,4,6};
	int first[288];
	int second[288];
	for(int k=0; k<4; k++)
	{
		int bpsc = bpsc_tbl[k];		// bpsc = number of coded bits per sub carrier
		int cbps = bpsc*48;			// cbps = number of coded bits per symbol
		int s = max(bpsc / 2, 1);
		for(int j = 0; j < cbps; j++) {
			first[j] = s * (j / s) + ((j + int(floor(16.0 * j / cbps))) % s);
		}

		for(int i = 0; i < cbps; i++) {
			second[i] = 16 * i - (cbps - 1) * int(floor(16.0 * i / cbps));
		}
		for(int i=0; i<cbps; i++)
			interleave_pattern[k][second[first[i]]] = i;
// 		for(int i=0; i<cbps; i++)
// 		{
// 			if (i%8 == 0)
// 				printf("\n");
// 			printf("%d,", interleave_pattern[k][i]);
// 		}
// 		printf("\n");
	}

	return 0;
}

uint32_t gettime()
{
	#ifdef WIN32
	return GetTickCount();
	#else
	struct timespec tv;  
	clock_gettime(CLOCK_MONOTONIC, &tv);    
	return (uint32_t)tv.tv_sec * 1000 + tv.tv_nsec/1000000;
	#endif
}

int * modulation2inerleaver_pattern(modulation mod)
{
	if (mod == BPSK)
		return interleave_pattern[0];
	if (mod == QPSK)
		return interleave_pattern[1];
	if (mod == QAM16)
		return interleave_pattern[2];
	if (mod == QAM64)
		return interleave_pattern[3];
	return NULL;
}

void descramble (uint8_t *decoded_bits, uint8_t *out_bytes, int bit_count) {

	int state = 0;

	for(int i = 0; i < 7; i++) {
		if(decoded_bits[i]) {
			state |= 1 << (6 - i);
		}
	}
	memset(out_bytes, 0, bit_count/8);
	out_bytes[0] = state;

	int feedback;
	int bit;

	for(int i = 7; i < bit_count; i++) {
		feedback = ((!!(state & 64))) ^ (!!(state & 8));
		bit = feedback ^ (decoded_bits[i] & 0x1);
		out_bytes[i/8] |= bit << (i%8);
		state = ((state << 1) & 0x7e) | feedback;
	}
}

int fft_complex(complex *in, complex *out, bool ifft = false)
{
	real_t I[64];
	real_t Q[64];

	real_t IO[64];
	real_t QO[64];

	for(int i=0; i<64; i++)
	{
		I[i] = in[i].real;
		Q[i] = in[i].image;
	}

	fft_real_t(64, ifft, I, Q, IO, QO);

	for(int i=0; i<64; i++)
	{
		out[i].real = IO[i];
		out[i].image = QO[i];
	}

	return 0;
}


int fft_complex_N(complex *in, complex *out, int N, bool ifft = false)
{
	real_t *I = new real_t[N];
	real_t *Q = new real_t[N];

	real_t *IO = new real_t[N];
	real_t *QO = new real_t[N];

	for(int i=0; i<N; i++)
	{
		I[i] = in[i].real;
		Q[i] = in[i].image;
	}

	fft_real_t(N, ifft, I, Q, IO, QO);

	for(int i=0; i<N; i++)
	{
		out[i].real = IO[i];
		out[i].image = QO[i];
	}

	delete [] I;
	delete [] Q;
	delete [] IO;
	delete [] QO;

	return 0;
}

complex LT_time_space[128];			// LT = long training, 2 repetition.
complex LT_frequency_space[64];		// one OFDM symbol only
complex ST_time_space[160];			// LT = long training, 2 repetition.
complex ST_frequency_space[64];		// one OFDM symbol only
complex LTS[64];

int init_training_sequence()
{
	// long training sequence is a predefined OFDM symbol
	// generate it by iFFT.
	float LTS[64] = 
	{
		0,	// DC
		1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,	// positive half
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,									// guard, both positive and negative halves
		1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,		// negative half
	};


	for(int i=0; i<64; i++)
	{
		LT_frequency_space[i].real = LTS[i] * (32767);
		::LTS[i].real = LTS[i];
	}

	fft_complex(LT_frequency_space, LT_time_space, true);
	for(int i=0; i<64; i++)
		LT_time_space[i+64] = LT_time_space[i];

	FILE * f = fopen("LT.pcm", "wb");
	for(int i=0; i<64; i++)
	{
		int16_t i16 = LT_time_space[i].real;
		int16_t q16 = LT_time_space[i].image;
		fwrite(&i16, 1, 2, f);
		fwrite(&q16, 1, 2, f);
	}

	fclose(f);


	// short training sequence
	int scale = 32767*sqrt(13.0/6.0);

	ST_frequency_space[4].real = -1 * scale;
	ST_frequency_space[4].image = -1 * -scale;

	ST_frequency_space[8].real = -1 * scale;
	ST_frequency_space[8].image = -1 * -scale;

	ST_frequency_space[12].real = 1 * scale;
	ST_frequency_space[12].image = 1 * -scale;

	ST_frequency_space[16].real = 1 * scale;
	ST_frequency_space[16].image = 1 * -scale;

	ST_frequency_space[20].real = 1 * scale;
	ST_frequency_space[20].image = 1 * -scale;

	ST_frequency_space[24].real = 1 * scale;
	ST_frequency_space[24].image = 1 * -scale;

	ST_frequency_space[64-4].real = 1 * scale;
	ST_frequency_space[64-4].image = 1 * -scale;

	ST_frequency_space[64-8].real = -1 * scale;
	ST_frequency_space[64-8].image = -1 * -scale;

	ST_frequency_space[64-12].real = -1 * scale;
	ST_frequency_space[64-12].image = -1 * -scale;

	ST_frequency_space[64-16].real = 1 * scale;
	ST_frequency_space[64-16].image = 1 * -scale;

	ST_frequency_space[64-20].real = -1 * scale;
	ST_frequency_space[64-20].image = -1 * -scale;

	ST_frequency_space[64-24].real = 1 * scale;
	ST_frequency_space[64-24].image = 1 * -scale;


	fft_complex(ST_frequency_space, ST_time_space, true);
	for(int i=0; i<16; i++)
		ST_time_space[i+64] = ST_time_space[i];
	for(int i=0; i<80; i++)
		ST_time_space[i+80] = ST_time_space[i];

	f = fopen("ST.pcm", "wb");
	for(int i=0; i<160; i++)
	{
		int16_t i16 = ST_time_space[i].real;
		int16_t q16 = ST_time_space[i].image;
		fwrite(&i16, 1, 2, f);
		fwrite(&q16, 1, 2, f);
	}

	fclose(f);

	return 0;
}

float constrain(float v, float _min, float _max)
{
	if (v < _min)
		return _min;
	if (v > _max)
		return _max;
	return v;
}

// assume a and b ranges from -PI to PI, 
float phase_sub(float a, float b)
{
	float v1 = a-b;
	float v2 = a+2*PI-b;
	float v3 = a-2*PI-b;

	v1 = fabs(v1)>fabs(v2) ? v2 : v1;
	return fabs(v1)>fabs(v3) ? v3 : v1;
}

int scrambler[127];

int init_scrambler()
{
	uint8_t start = 0x7f;
	uint8_t a = start;
	for(int i = 1;; i++) {
// 		printf("%x", a>>7);
		scrambler[i-1]=a>>7;
// 		if (i%8==0)
// 			printf(" ");
		int newbit = (((a >> 1) ^ (a >> 4)) & 1);
		a = ((a >> 1) | (newbit<<7)) & 0xff;
		if (a == start) {
			printf("\n");
			break;
		}

	}

	return 0;
}

int depuncture(uint8_t *bits, int input_count, puncturing pun)
{
	if (pun == _1_2)
		return input_count;
	if (pun == _2_3)
	{
		int output_count = input_count * 4 / 3;
		int block_count = input_count/3;
		uint8_t *v = new uint8_t[input_count];
		memcpy(v, bits, input_count);

		for(int i=block_count-1; i>=0; i--)
		{
			bits[i*4+0] = v[i*3+0];
			bits[i*4+1] = v[i*3+1];
			bits[i*4+2] = v[i*3+2];
			bits[i*4+3] = 2;
		}

		delete v;
		return output_count;
	}

	if (pun == _3_4)
	{
		int block_count = input_count/4;
		int output_count = block_count * 6;
		uint8_t *v = new uint8_t[input_count];
		memcpy(v, bits, input_count);

		for(int i=block_count-1; i>=0; i--)
		{
			bits[i*6+0] = v[i*4+0];
			bits[i*6+1] = v[i*4+1];
			bits[i*6+2] = v[i*4+2];
			bits[i*6+3] = 2;
			bits[i*6+4] = 2;
			bits[i*6+5] = v[i*4+3];
		}

		delete v;

		return output_count;
	}

	return 0;
}

int puncture(uint8_t *bits, int input_count, puncturing pun)
{
	if (pun == _1_2)
		return input_count;
	if (pun == _2_3)
	{
		int block_count = input_count/4;
		int output_count = block_count*3;

		for(int i=0; i<block_count; i++)
		{
			bits[i*3+0] = bits[i*4+0];
			bits[i*3+1] = bits[i*4+1];
			bits[i*3+2] = bits[i*4+2];
		}

		return output_count;
	}

	if (pun == _3_4)
	{
		int block_count = input_count / 6;
		int output_count = block_count * 4;

		for(int i=0; i<block_count; i++)
		{
			bits[i*4+0] = bits[i*6+0];
			bits[i*4+1] = bits[i*6+1];
			bits[i*4+2] = bits[i*6+2];
			bits[i*4+3] = bits[i*6+5];
		}

		return output_count;
	}

	return 0;
}

int ones[256];
void init_ones()
{
	for(int n=0; n<256; n++)
	{
		ones[n] = 0;
		for(int i = 0; i < 8; i++) {
			if(n & (1 << i)) {
				ones[n]++;
			}
		}
	}
}

void convolutional_encoding(const uint8_t *in, uint8_t *out, int input_count)
{
	uint8_t state = 0;

	for(int i = 0; i < input_count; i++) 
	{
		assert(in[i] == 0 || in[i] == 1);
		state = ((state&0x3f) << 1) | in[i];			// shift left and input data on LSB
		out[i * 2 + 0] = ones[state & 0x6d] & 1;		// 802.11 generator polynomial A: binary 1101101(0x6d, 0155)
		out[i * 2 + 1] = ones[state & 0x4f] & 1;		// 802.11 generator polynomial B: binary 1001111(0x4f, 0117)
	}
}

int test_convolutional_code()
{
	init_ones();

	uint8_t data[] = {1, 0, 1, 0, 1, 1, 1, 1, 1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0};

	char source[512] = {0};
	uint8_t encoded[512] = {0};
	uint8_t decoded[512] = {0};

	printf("data:\n");
	for(int i=0; i<sizeof(data); i++)
		printf("%d ", data[i]);
	printf("\n");

	convolutional_encoding(data, encoded, sizeof(data));

	printf("encoded:\n");
	for(int i=0; i<sizeof(data); i++)
		printf("%d%d ", encoded[i*2+0], encoded[i*2+1]);
	printf("\n");


	puncturing pun = _3_4;
	int s1 = puncture(encoded, sizeof(data)*2, pun);
// 	for(int i=0; i<14; i++)
// 		encoded[i] ^= 1;
// 
// 	printf("corrupted:\n");
// 	for(int i=0; i<sizeof(data); i++)
// 		printf("%d%d ", encoded[i*2+0], encoded[i*2+1]);
// 	printf("\n");
	int s2 = depuncture(encoded, s1, pun);

	assert(s2 == sizeof(data)*2);

	viterbi_decoder dec;
	dec.decode((uint8_t*)encoded, decoded, (sizeof(data)+15)/16*16, 9);

	printf("decoded:\n");
	for(int i=0; i<sizeof(data); i++)
		printf("%d%s", decoded[i], decoded[i] == data[i] ? " " : "x");
	printf("\n");

	return 0;
}

int frame_decoding(complex * s, int sample_count, float frequency_offset, uint8_t *out_data, int *valid_data_len)
{
	int l = gettime();
	*valid_data_len = 0;
	FILE * f;

	// apply frequency offset for preambles and "SIGNAL" symbol.
 	int MAX_SIGNAL_FIELD_POS = 160+160+80;
	complex_frequency_offset(s, MAX_SIGNAL_FIELD_POS, frequency_offset);

	// OFDM symbol alignment
	// use long training sequence to find first symbol
// 	f =fopen("alignment.csv", "wb");
	//fprintf(f, "N,P,A\n");
	float correction_value[3] = {0};
	int correction_pos[3] = {0};
	for(int n=0; n<MAX_SIGNAL_FIELD_POS-80-64; n++)
	{
		float corelation = correalation((float*)&s[n], (float*)LT_time_space, 128);

//		fprintf(f, "%d,%f,%f\n", n, corelation, s[n].magnitude());
		for(int i=0; i<3; i++)
		{
			if (corelation > correction_value[i])
			{
				for(int j=i+1; j<3; j++)
				{
					correction_value[j] = correction_value[j-1];
					correction_pos[j] = correction_pos[j-1];
				}

				correction_value[i] = corelation;
				correction_pos[i] = n;

				break;
			}
		}
	}

	int peak_pos = 0;
	for(int i=0; i<3; i++)
		if (correction_pos[i] > peak_pos)
			peak_pos = correction_pos[i];
// 	fclose(f);

	printf("long training sequence @ %d of preamble\n", peak_pos);
	int symbol_start = peak_pos + 64;
	printf("symbols start @ %d (%.1fus)\n", symbol_start, symbol_start/20.0f);

	// fine frequency offset
	complex _df;
	for(int i=0; i<64; i++)
		_df += s[symbol_start-64+i] * s[symbol_start-128+i].conjugate();
	float df = _df.argument()/64;

	printf("fine frequency offset:%f degree / sample (%f Mhz)\n", df * 180 / PI, 20 / (2 * PI / df) );

	// apply fine frequency offset of header and add to coarse frequency for main body
	complex_frequency_offset(s, MAX_SIGNAL_FIELD_POS, -df);
	frequency_offset -= df;


	// initial phase and amplitude equalization
	static complex lt_rx_fft[64];
	fft_complex(s+symbol_start-64, lt_rx_fft);

	static complex h[64];
	for(int i=-26; i<=26; i++)
	{
		int n = i > 0 ? i : (i+64);

		h[n] = LTS[n] / lt_rx_fft[n];
	}


	// decode "SIGNAL" symbol, BPSK, 1/2 convolutional coded
	// if this symbol fails to decode, we drop the whole frame
	static complex signal_fft[64];
	fft_complex(&s[symbol_start+16], signal_fft);
	ALIGN uint8_t signal_bits[640] = {0};
	ALIGN uint8_t signal_bits_deinterleaved[640] = {0};
	ALIGN uint8_t signal_decoded_bits[640];
	int signal_pos = 0;
	printf("signal bits:\n");
	for(int i=-26; i<=26; i++)
	{
		if (i == 0 || i == 7 || i == 21 || i == -7 || i== -21)
			continue;
		int n = i>0?i:(i+64);
		signal_fft[n] = signal_fft[n] * h[n];
		int bit = signal_fft[n].real > 0 ? 1 : 0;
		signal_bits[signal_pos++] = bit;

		printf("%d", bit);
		if (signal_pos %8 == 0)
			printf(" ");
	}

	for(int i=0; i<48; i++)
		signal_bits_deinterleaved[i] = signal_bits[interleave_pattern[0][i]];

	printf("\ndecoded:\n");
	viterbi_decoder dec;
	dec.decode(signal_bits_deinterleaved, signal_decoded_bits, 48);
	for(int i=0; i<24; i++)
	{
		if (i % 8 == 0 && i > 0)
			printf(" ");
		printf("%d", signal_decoded_bits[i]);
	}
	printf("\n");

	int rate_code = (signal_decoded_bits[0] << 3) | (signal_decoded_bits[1] << 2) | (signal_decoded_bits[2] << 1) | signal_decoded_bits[3];

	printf("rate:%dMbps\n", rate_code_to_mbps(rate_code));
	printf("reserved bytes:%d\n", signal_decoded_bits[4]);
	int length = 0;
	for(int i=0; i<12; i++)
		length |= signal_decoded_bits[i+5] << i;
	printf("length:%d\n", length);
	int parity = 0;
	for(int i=0; i<17; i++)
		parity ^= signal_decoded_bits[i];
	printf("parity:%d%s%d\n", parity, parity == signal_decoded_bits[17] ? "=" : "!=", signal_decoded_bits[17]);
	printf("tail:");
	bool tail_is_zero = true;
	for(int i=0; i<6; i++)
	{
		printf("%d", signal_decoded_bits[i+18]);
		tail_is_zero &= !signal_decoded_bits[i+18];
	}
	printf("%s\n", tail_is_zero?"(OK)":"(ERROR)");

	if (rate_code_to_mbps(rate_code) == 0
		|| signal_decoded_bits[4] != 0
		|| parity != signal_decoded_bits[17]
		|| !tail_is_zero
		)
	{
		printf("invalid signal field\n");
		return symbol_start+80;
	}

	int rate = rate_code_to_mbps(rate_code);
	modulation mod = rate2modulation(rate);
	puncturing pun = rate2puncturing(rate);
	printf("modulation:%s\n", modulation_name(mod));
	printf("coding rate:%s\n", puncturing_name(pun));

	int data_bits_per_symbol = mod*48 * pun/12;

	int data_symbol_count = (length*8+ 16+6 + data_bits_per_symbol-1) / (data_bits_per_symbol);		// *8: 8bits per byte, +16: service field, +6: tailing
	printf("%d data symbols\n", data_symbol_count);

	int symbol_count = data_symbol_count + 1;		// 1 = signal symbol

	if (sample_count - symbol_start < symbol_count * 80)
	{
		printf("not enough data sample");
		return sample_count;
	}

	// apply frequency offset for data symbols
	complex_frequency_offset(s+MAX_SIGNAL_FIELD_POS,
		symbol_start + symbol_count * 80-MAX_SIGNAL_FIELD_POS, frequency_offset, MAX_SIGNAL_FIELD_POS*frequency_offset);


	// fft and qualize all symbols
	static complex symbols[200][64];

	for(int i=0; i<symbol_count; i++)
	{
		fft_complex(s+symbol_start+i*80+16, symbols[i]);

		// equalizing
		for(int j=0; j<64; j++)
			symbols[i][j] = symbols[i][j] * h[j];
	}


	// compensate for sampling rate error, which introduce a triangular shaped phase shift to all subcarriers
	// use pilot tones to compensate for it
	float pilot_phase[4][200];
	float pilot_phase_error[4][200];
	for(int i=0; i<symbol_count; i++)
	{
		pilot_phase[0][i] = symbols[i][-21+64].argument();
		pilot_phase[1][i] = symbols[i][-7+64].argument();
		pilot_phase[2][i] = symbols[i][7].argument();
		pilot_phase[3][i] = symbols[i][21].argument();

		for(int j=0; j<4; j++)
		{
			int bit = j==3 ? 1-scrambler[i%127] : scrambler[i%127];
			float bit_phase = bit ? PI : 0;
			pilot_phase_error[j][i] = phase_sub(pilot_phase[j][i], bit_phase);
		}
	}

	float k1[4], b1[4];
	for(int i=0; i<4; i++)
	{
		k1[i] = 0;
		for(int j=1; j<symbol_count; j++)
		{
			float diff = phase_sub(pilot_phase_error[i][j], pilot_phase_error[i][j-1]);
// 			printf("%f\n", diff);
			k1[i] += diff;
		}

		k1[i] /= symbol_count-1;
		b1[i] = pilot_phase_error[i][0];
	}

	float phase_error[200][64];
	float k[4],b[4];
	for(int i=0; i<4; i++)
	{
		line_fitting_y(pilot_phase_error[i], symbol_count, k+i, b+i);

// 		printf("%f*i+%f\t%f*i+%f\n", k[i], b[i], k1[i], b1[i]);

		for(int j=0; j<symbol_count; j++)
			pilot_phase_error[i][j] = k[i]*j+b[i];
	}

	float idx[4] = {-21, -7, 7, 21};
	float ki,bi;
	line_fitting(idx, k1, 4, &ki,&bi);

	for(int j=0; j<symbol_count; j++)
	{
		for(int i=-26; i<=26; i++)
		{
			if (i == 0)
				continue;

			int idx = i>0 ? i : i+64;
			float k_i = ki*i+bi;

// 			printf("%d,%f\n", i, k_i);
			
			float phase_error = -(k_i * j + (b1[0]+b1[1]+b1[2]+b1[3])/4);
			symbols[j][idx] = symbols[j][idx] * cordic(phase_error);//offset;
		}
	}

	// now decode the data symbols
	ALIGN uint8_t service_and_data_bits[8192*8] = {0};
	ALIGN uint8_t service_and_data_deinterleaved_bits[8192*8] = {0};
	ALIGN uint8_t service_and_data_decoded_bits[8192*8];
	int pos = 0;

	mapper demapper = modulation2demapper(mod);
	mapper remapper = modulation2mapper(mod);
	static complex remapped_symbol[64];
	int * pattern = modulation2inerleaver_pattern(mod);
	int bits_per_symbol = 48*mod;
	int data_bit_per_symbol = bits_per_symbol * pun / 12;

	float EVM = 0;
	float sum = 0;
	float sum2 = 0;

	static complex h2[64];
	for(int i=0; i<64; i++)
	{
		h[i].real = 1;
		h[i].image = 0;
	}

	float alpha = 0.3;

	for(int i=0; i<data_symbol_count; i++)
	{
		// continual equalizing
		for(int j=-26; j<26; j++)
		{
			if (j == 0 || j == -21 || j == 21 || j == 7 || j == -7)	// ignore DC/pilots
				continue;
			int n = j<0 ? j+64 : j;

			symbols[i+1][n] *= h[n];
		}

		// map bits
		uint8_t *p = service_and_data_bits + bits_per_symbol*i;
		demapper(symbols[i+1], p);
		remapper(remapped_symbol, p);

		// EVM and h calculation
		for(int j=-26; j<26; j++)
		{
			int n = j<0 ? j+64 : j;
			if (j == 0)										// ignore DC
				continue;

			if (j == -21 || j == 21 || j == 7 || j == -7)	// calculate SNR of pilots
			{
				sum += symbols[i+1][n].magnitude();
				sum2 += symbols[i+1][n].sq_magnitude();
				continue;
			}

			complex _h = remapped_symbol[n] / symbols[i+1][n];
			h[n] = h[n]*(1-alpha) + _h * alpha;
			float _EVM = (symbols[i+1][n] - remapped_symbol[n]).magnitude();
			EVM += _EVM;
		}

		// deinterleave
		uint8_t *p2 = service_and_data_deinterleaved_bits+ bits_per_symbol*i;
		for(int j=0; j<bits_per_symbol; j++)
			p2[j] = p[pattern[j]];
	}


// 	FILE * constellation = fopen("constellation.csv", "wb");
// 	fprintf(constellation, "N,I,Q,P,A\n");
// 	for(int i=1; i<symbol_count; i++)
// 	{
// 		for(int j=-26; j<=26; j++)
// 		{
// 			if (j == 0)
// 				continue;
// 
// 			int n = j>0?j:j+64;
// 			if (!(n==7||n==21||n==-7||n==-21||n==64-7||n==64-21))
// 				continue;
// 
// 			complex v = symbols[i][n];
// 
// 			fprintf(constellation, "%d,%f,%f,%f,%f\n", i, v.real, v.image, v.argument(), v.magnitude());
// 		}
// 	}
// 	fclose(constellation);

	EVM /= 52*data_symbol_count;

	float avg = sum / (4*data_symbol_count);
	float variance = sum2 / (4*data_symbol_count) - avg*avg;
	float SNR = avg*avg/variance;
	float SNRdb = 10*log10(SNR);
	printf("EVM>=%.1f%%(%.1fdb)\n", EVM*100, log10(EVM)*20);
	printf("SNR<=%.1fdb\n", SNRdb);

	printf("sample_count = %d\n", sample_count);
	fflush(stdout);
	int bits_count = depuncture(service_and_data_deinterleaved_bits, bits_per_symbol*data_symbol_count, pun);
	printf("bits_count = %d\n", bits_count);
	fflush(stdout);
	int ntraceback = 5;
	if (pun == _2_3)
		ntraceback = 9;
	if (pun == _3_4)
		ntraceback = 10;
	dec.decode(service_and_data_deinterleaved_bits, service_and_data_decoded_bits, bits_count+16, ntraceback);
	uint8_t out_bytes[8192];
	descramble(service_and_data_decoded_bits, out_bytes, data_bit_per_symbol*data_symbol_count);

	uint32_t crc = crc32_80211(out_bytes+2, length-4);
	uint8_t *p = out_bytes+ 2+length-4;
	uint32_t crc_received = (p[3] << 24) | (p[2] << 16) | (p[1] << 8) | p[0];
	printf("FCS(calculated)=0x%08x\n", crc);
	printf("FCS(received)=0x%08x\n", crc_received);

// 	f = fopen("data.bin", "wb");
// 	fwrite(out_bytes+2, 1, length, f);
// 	fclose(f);

	if (crc == crc_received)
	{
		memcpy(out_data, out_bytes+2, length);
		*valid_data_len =  length;
	}
	else
	{
		*valid_data_len = 0;
	}

	printf("%dms\n", gettime()-l);
	return symbol_start + symbol_count * 80;
}


int tx(uint8_t *psdu, int count, complex **out, int mbps = 6)
{
	float scale = 32767;

	modulation mod = rate2modulation(mbps);
	puncturing pun = rate2puncturing(mbps);
	int data_bits_per_symbol = mod*48 * pun/12;
	int bits_per_symbol = 48*mod;
	int data_symbol_count = (count*8+16+6 + data_bits_per_symbol-1) / data_bits_per_symbol;		// 16: SERVICE field, 6: tail
	int data_bits_count_padded = data_symbol_count * data_bits_per_symbol;
	int sample_count = (data_symbol_count+1)*80+320+200;
	complex *s = new complex[sample_count];
	*out = s;

	// short training sequence
	for(int i=0; i<160; i++)
		s[i] = ST_time_space[i];

	// GI2
	for(int i=0; i<32; i++)
		s[i+160] = LT_time_space[i+32];

	// long training sequence
	for(int i=0; i<128; i++)
		s[i+192] = LT_time_space[i];

	int rate_code = mbps_to_rate_code(mbps);
	uint8_t signal_bits[24] = 
	{
		(rate_code>>3)&1, (rate_code>>2)&1, (rate_code>>1)&1, (rate_code>>0)&1,
		0,				// "reserved bit"
	};

	for(int i=0; i<12; i++)						// "LENGTH"
		signal_bits[i+5] = (count >> i) & 1;

	for(int i=0; i<17; i++)						// "parity"
		signal_bits[17] ^= signal_bits[i];

	// encode ( no puncturing for signal symbol)
	uint8_t signal_bits_encoded[48];
	convolutional_encoding(signal_bits, signal_bits_encoded, 24);

	// interleave
	uint8_t signal_bits_interleaved[48];
	for(int i=0; i<48; i++)
		signal_bits_interleaved[interleave_pattern[0][i]] = signal_bits_encoded[i];

	// modulation
	complex signal_symbol[64];
	for(int i=-26, j=0; i<=26; i++)
	{
		if (i == 0 || i == -7 || i == -21 || i == 7 || i == 21)
			continue;
		
		signal_symbol[i>0?i:i+64] = signal_bits_interleaved[j++] ? scale : - scale;
	}

	// pilot tones
	signal_symbol[7].real = scale;
	signal_symbol[21].real = -scale;
	signal_symbol[64-7].real = scale;
	signal_symbol[64-21].real = scale;


	complex signal_symbol_ifft[64];
	fft_complex(signal_symbol, signal_symbol_ifft, true);



	for(int i=0; i<16; i++)
		s[i+320] = signal_symbol_ifft[i+48];
	for(int i=0; i<64; i++)
		s[i+336] = signal_symbol_ifft[i];

	// data scrambling and add service field
	uint8_t data_bits[4096*8] = {1, 1, 1, 1, 1, 1, 1, 0};		// we always use scrambler starting from all ones state
	for(int i=0; i<9; i++)
	{
		data_bits[i+7] ^= scrambler[i];
	}


	for(int i=0; i<count; i++)
	{
		for(int j=0; j<8; j++)
		{
			data_bits[16+i*8+j] = ((psdu[i]>>j)&1) ^ (scrambler[(i*8+j+9)%127]);
		}
	}

// 	uint8_t descramble_bytes[4096] = {0};
// 	descramble(data_bits, descramble_bytes, data_bits_count_padded);

	// convolutional encoding and puncturing
	uint8_t data_bits_encoded[8192*8];
	convolutional_encoding(data_bits, data_bits_encoded, data_bits_count_padded);
	int bits_count_punctured = puncture(data_bits_encoded, data_bits_count_padded*2, pun);
	assert(bits_count_punctured == data_symbol_count*bits_per_symbol);

	// interleave
	uint8_t data_bits_interleaved[8192*8];
	int *pattern = modulation2inerleaver_pattern(mod);

	for(int i=0; i<data_symbol_count; i++)
	{
		uint8_t *p2 = data_bits_interleaved+ bits_per_symbol*i;
		uint8_t *p = data_bits_encoded+ bits_per_symbol*i;
		for(int j=0; j<bits_per_symbol; j++)
			p2[pattern[j]] = p[j];
	}

	// modulation, add pilot tones and output
	mapper mapper = modulation2mapper(mod);
	for(int i=0; i<data_symbol_count; i++)
	{
		complex symbol_fre[64];
		complex symbol_time[64];

		// data mapper
		mapper(symbol_fre, data_bits_interleaved + bits_per_symbol*i);
		
		// pilot tone
		float pilot = scrambler[(i+1)%127] ? -1 : 1;
		symbol_fre[7].real = pilot;
		symbol_fre[21].real = -pilot;
		symbol_fre[64-7].real = pilot;
		symbol_fre[64-21].real = pilot;

		// scale
		for(int j=0; j<64; j++)
		{
			symbol_fre[j].real *= scale;
			symbol_fre[j].image *= scale;
		}

		// ifft
		fft_complex(symbol_fre, symbol_time, true);

		// out: GI + symbol
		for(int j=0; j<16; j++)
			s[400+i*80+j] = symbol_time[j+48];		// 400: ST+LT+SIGNAL
		for(int j=0; j<64; j++)
			s[400+i*80+16+j] = symbol_time[j];		// 400: ST+LT+SIGNAL
	}


	FILE * f = fopen("out_8bit.pcm", "wb");
	for(int i=0; i<sample_count; i++)
	{
		int8_t I = -s[i].real / 256;
		int8_t Q = -s[i].image / 256;
		fwrite(&Q, 1, 1, f);
		fwrite(&I, 1, 1, f);
	}
	fclose(f);

	f = fopen("out_16bit.pcm", "wb");
	for(int i=0; i<sample_count; i++)
	{
		int m = 1<<6;
		int16_t I = int((s[i].real + (rand()&0xff-127)*m/float(m*2)) / m) * m;
		int16_t Q = int((s[i].image + (rand()&0xff-127)*m/float(m*2)) / m) * m;
		fwrite(&Q, 1, 2, f);
		fwrite(&I, 1, 2, f);

	}
	fclose(f);

// 	uint8_t decoded_data[4096] = {0};
// 	int decoded_count = 0;
// 	frame_decoding(s+80, sample_count, 0.001, decoded_data, &decoded_count);

	return sample_count;
}





int16_t queue[256*1024 + MAX_SYMBOLS_PER_SAMPLE * 80*2 + 500] = {0};
int16_t a[256*1024 + MAX_SYMBOLS_PER_SAMPLE * 80*2 + 500] = {0};
int16_t p[128*1024 + MAX_SYMBOLS_PER_SAMPLE * 80*2 + 500] = {0};
int queue_count = 0;

int slice(int16_t *new_data, int count, bool flush = false)
{
	// overflow check
	if (count > 256*1024)
		return -1;

	// now padding is enough
	if (new_data && count)
	{
		memcpy(queue+queue_count, new_data, count*2);		
		queue_count += count;
	}

	// SIMD complex multiply / amplitude
	complex16_abj_amp(queue, queue+32, a, p, queue_count-32);

	static FILE * f = fopen("Z:\\auto.pcm", "wb");
	int Nwindow = 48;

	int moving_avg_a[2] = {0};
	int moving_avg_p = 0;
	for(int i=0; i<Nwindow; i++)
	{
		moving_avg_p += p[i];
		moving_avg_a[0] += a[i*2];
		moving_avg_a[1] += a[i*2+1];
	}

	int counter = 0;

	static FILE * fout = fopen("data.pkt", "wb");
	for(int i=Nwindow; i<(flush?(queue_count+MAX_SYMBOLS_PER_SAMPLE * 80*2 + 500):queue_count)/2-MAX_SYMBOLS_PER_SAMPLE*80; i++)
	{
		int auto_denum = (moving_avg_p/100*moving_avg_p/640);
		if (auto_denum > 0)
		{
			int auto10000 = (moving_avg_a[0]/64*moving_avg_a[0]+moving_avg_a[1]/64*moving_avg_a[1]) / auto_denum;
			static int threshold = (0.8*0.8)*1000;

			if (auto10000 > threshold)
				counter ++;
			else
				counter = 0;

			if (counter == 48)
			{
				int complex_offset[2] = {moving_avg_a[0], moving_avg_a[1]};
				for(int j=0; j<48; j++)
				{
					complex_offset[0] += a[(i+j)*2+0];
					complex_offset[1] += a[(i+j)*2+1];
				}

				float offset = atan2((double)complex_offset[1], (double)complex_offset[0])/16;
				printf("coarse frequency offset: %f degress/sample, %fMhz\n", offset*180/PI, 20*offset / (2 * PI));

				// copy to complex and apply coarse frequency offset;
				static complex pkt[MAX_SYMBOLS_PER_SAMPLE*80];
				int16_t pkt2[MAX_SYMBOLS_PER_SAMPLE*80*2];
				int pkt_size = MAX_SYMBOLS_PER_SAMPLE*80;

// 				complex16_frequency_offset(queue+i*2, pkt2, pkt_size*2, offset);
// 
// 				offset = 0;

				complex_from_int16(pkt, queue+i*2, pkt_size, true);

				int valid_data_len = 0;
				uint8_t out_data[4096];



				int pkt_result = frame_decoding(pkt, pkt_size, -offset, out_data, &valid_data_len);
				

				static int pkt_count = 0;
				static int valid_pkt_count = 0;
				pkt_count++;
				if (valid_data_len > 0)
				{
					fwrite(&valid_data_len, 1, 4, fout);
					fwrite(out_data, 1, valid_data_len, fout);

					valid_pkt_count++;
				}

				fprintf(stderr, "\r%d/%d pkt", valid_pkt_count, pkt_count);


			}

		}

		moving_avg_p -= p[i-Nwindow];
		moving_avg_a[0] -= a[(i-Nwindow)*2];
		moving_avg_a[1] -= a[(i-Nwindow)*2+1];

		moving_avg_p += p[i];
		moving_avg_a[0] += a[i*2];
		moving_avg_a[1] += a[i*2+1];
	}

	if (flush)
		fflush(fout);

	if (queue_count > (MAX_SYMBOLS_PER_SAMPLE*80)*2)
	{
		int dn = queue_count-MAX_SYMBOLS_PER_SAMPLE*80*2;

		memmove(queue, queue+dn, MAX_SYMBOLS_PER_SAMPLE*80*4);

		assert(queue_count-dn == MAX_SYMBOLS_PER_SAMPLE*80*2);
		queue_count = MAX_SYMBOLS_PER_SAMPLE*80*2;
	}

	return count;
}

int main()
{
	init_training_sequence();
	init_scrambler();
	init_interleaver_pattern();
	init_cordic();
	init_ones();

	uint8_t psdu[500] = 
	{
		0x08, 0x01, 						// FC: frame control
		0x00, 0x01,							// DID: duration or ID
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,	// address1
		0x13, 0x22, 0x33, 0x44, 0x55, 0x66,	// address2
		0x13, 0x22, 0x33, 0x44, 0x55, 0x66,	// address3
		0x10, 0x86,							// SC: Sequence control
	};


	uint32_t crc = crc32_80211(psdu, sizeof(psdu)-4);
	psdu[sizeof(psdu)-4+3] = (crc >> 24)&0xff;
	psdu[sizeof(psdu)-4+2] = (crc >> 16)&0xff;
	psdu[sizeof(psdu)-4+1] = (crc >> 8)&0xff;
	psdu[sizeof(psdu)-4+0] = (crc >> 0)&0xff;
	complex *s_gen = NULL;
	tx(psdu, sizeof(psdu), &s_gen, 48);
	delete [] s_gen;

	char file[] = "testcases/wifi_rt5572_usrpb205_18mbps_nonlinear.pcm";
	bool is8bit = strstr(file, "8bit");
	FILE * f = fopen(file, "rb");
	if (!f)
	{
		printf("failed opening file\n");
		return 0;
	}

	int16_t *data = new int16_t[128*1024];
	while(!feof(f))
	{
		int byte_per_sample = is8bit?1:2;
		int count = fread(data, 1, 128*1024*byte_per_sample, f)/byte_per_sample;
		if (is8bit)
		{
			int8_t *p = (int8_t*)data;
			for(int i=count-1; i>=0; i--)
				data[i] = p[i]*256;
		}

		slice(data, count);
	}

	delete []data;
	fclose(f);
	slice(NULL, 0, true);

	return 0;
}
