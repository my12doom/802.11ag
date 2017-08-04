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

#ifdef WIN32
#include <Windows.h>
#endif

using namespace gr::ieee802_11;

#define countof(x) (sizeof(x)/sizeof(x[0]))

const int CORDIC_TBL_COUNT = 100000;
const int MAX_SYMBOLS_PER_SAMPLE = 200;
complex cordic_tbl[CORDIC_TBL_COUNT];
int interleave_pattern[4][288];	// [BPSK, QPSK, 16QAM, 64QAM] [max: 6*48bits]

int init_cordic()
{
	for(int i=0; i<CORDIC_TBL_COUNT; i++)
	{
		float phase = i*2*PI/CORDIC_TBL_COUNT;
		cordic_tbl[i] = complex(cos(phase), sin(phase));
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

uint8_t map_BPSK(complex in)
{
	return in.real > 0 ? 1 : 0;
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
		LT_frequency_space[i].real = LTS[i] * (32767);

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

int frame_decoding(complex * s, int sample_count, uint8_t *out_data, int *valid_data_len)
{
	int l = gettime();
	*valid_data_len = 0;

	FILE * f;

	// OFDM symbol alignment
	// use long training sequence to find first symbol
	f =fopen("alignment.csv", "wb");
	fprintf(f, "N,P,A\n");
	float autocorrection_value[3] = {0};
	int autocorrection_pos[3] = {0};
	for(int n=0; n<360+50; n++)			// 160: number of sample of short training sequence
		// 320: number of sample of long training sequence
	{
		complex a;
		float p = 0;
		for(int k=0; k<64; k++)
		{
			a += s[n+k] * LT_time_space[k].conjugate();
			p += s[n+k].sq_magnitude();//(s[preamble_start+n+k] * s[preamble_start+n+k].conjugate()).real;
		}

		p = max(max(p, s[n+64].magnitude()), s[n].magnitude());

		float autocorelation = a.magnitude()/p;
		fprintf(f, "%d,%f,%f\n", n, autocorelation, s[n].magnitude());

		for(int i=0; i<3; i++)
		{
			if (autocorelation > autocorrection_value[i])
			{
				for(int j=i+1; j<3; j++)
				{
					autocorrection_value[j] = autocorrection_value[j-1];
					autocorrection_pos[j] = autocorrection_pos[j-1];
				}

				autocorrection_value[i] = autocorelation;
				autocorrection_pos[i] = n;

				break;
			}
		}
	}

	int peak_pos = 0;
	for(int i=0; i<3; i++)
		if (autocorrection_pos[i] > peak_pos)
			peak_pos = autocorrection_pos[i];
	fclose(f);

	printf("long training sequence @ %d of preamble\n", peak_pos);
	int symbol_start = peak_pos + 64;
	printf("symbols start @ %d (%.1fus)\n", symbol_start, symbol_start/20.0f);

	// fine frequency offset
	complex _df;
	for(int i=0; i<64; i++)
	{
		complex _f = s[symbol_start-64+i] * s[symbol_start-128+i].conjugate();
		_df += _f;

	}

	float df = _df.argument()/64;

	printf("fine frequency offset:%f degree / sample (%f Mhz)\n", df * 180 / PI, 20 / (2 * PI / df) );


	// initial phase and amplitude equalization
	complex lt_rx_fft[64];
	fft_complex(s+symbol_start-64, lt_rx_fft);
// 	fft_complex(LT, lt_fft);

	complex h[64];
	for(int i=-26; i<=26; i++)
	{
		if (i == 0)
			continue;

		int n = i > 0 ? i : (i+64);
		//printf("LT[%d]=%.1f,%.1f, %.2f\n", i, lt_fft[i].real, lt_fft[i].image, lt_fft[i].argument());

		float phase_diffs = phase_sub(LT_frequency_space[n].argument(), lt_rx_fft[n].argument());
		float amplitude_normalizer = 1/lt_rx_fft[n].magnitude();

		h[n] = complex::from_phase_magnitude(phase_diffs, amplitude_normalizer);

// 		assert(abs(phase_diffs- (LT_frequency_space[n]*lt_rx_fft[n].conjugate()).argument()) < 0.05);
		complex normalized = lt_rx_fft[n]*h[n];

// 		printf("%d,%.2f,%.4f\n", i, normalized.argument(), normalized.magnitude());
	}


	// "SIGNAL" symbol, BPSK, 1/2 convolutional coded
	// if this symbol fails to decode, we drop the whole frame
	complex signal_fft[64];
	fft_complex(&s[symbol_start+16], signal_fft);
	ALIGN uint8_t signal_bits[64] = {0};
	ALIGN uint8_t signal_bits_deinterleaved[64] = {0};
	ALIGN uint8_t signal_decoded_bits[64];
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
	viterbi_decoder *dec = new viterbi_decoder;
	dec->decode(signal_bits_deinterleaved, signal_decoded_bits, 48);
	delete dec;
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
// 		|| length > 1024
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


	// fft all symbols
	complex symbols[200][64];

	for(int i=0; i<symbol_count; i++)
	{
		fft_complex(s+symbol_start+i*80+16, symbols[i]);

		// equalizing
		float phase_offset_per_symbol = -df * 80;
		float phase_offset = (i+1)*phase_offset_per_symbol;
// 		complex offset(cos(phase_offset), sin(phase_offset));
// 		complex offset2 = cordic(phase_offset);
		for(int j=0; j<64; j++)
			symbols[i][j] = symbols[i][j] * h[j] * cordic(phase_offset);//offset;
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

	for(int i=-26; i<=26; i++)
	{
		if (i == 0)
			continue;

		int idx = i>0 ? i : i+64;
		float k_i = ki*i+bi;

// 		printf("%d,%f\n", i, k_i);
		
		for(int j=0; j<symbol_count; j++)
		{
			float phase_error = -(k_i * j + (b1[0]+b1[1]+b1[2]+b1[3])/4);
// 			complex offset(cos(phase_error), sin(phase_error));
			symbols[j][idx] = symbols[j][idx] * cordic(phase_error);//offset;
		}
	}

	// now the data symbols
	ALIGN uint8_t service_and_data_bits[8192*8] = {0};
	ALIGN uint8_t service_and_data_deinterleaved_bits[8192*8] = {0};
	ALIGN uint8_t service_and_data_decoded_bits[8192*8];
	int pos = 0;

	FILE * constellation = fopen("constellation.csv", "wb");
	fprintf(f, "N,P,A\n");
	for(int i=1; i<symbol_count; i++)
	{
		for(int j=-26; j<=26; j++)
		{
			if (j == 0)
				continue;

			int n = j>0?j:j+64;
			if ((n==7||n==21||n==-7||n==-21||n==64-7||n==64-21))
				continue;

			complex v = symbols[i][n];

			fprintf(constellation, "%d,%f,%f,%f\n", i, v.real, v.image, v.argument());
		}
	}
	fclose(constellation);

	mapper demapper = modulation2demapper(mod);
	int * pattern = modulation2inerleaver_pattern(mod);
	int bits_per_symbol = 48*mod;
	int data_bit_per_symbol = bits_per_symbol * pun / 12;

	for(int i=0; i<data_symbol_count; i++)
	{
		// map bits
		uint8_t *p = service_and_data_bits + bits_per_symbol*i;
		demapper(symbols[i+1], p);

		// deinterleave
		uint8_t *p2 = service_and_data_deinterleaved_bits+ bits_per_symbol*i;
		for(int j=0; j<bits_per_symbol; j++)
			p2[j] = p[pattern[j]];
	}

// 	for(int i=0; i<bits_per_symbol*data_symbol_count; i++)
// 	{
// 		printf("%d", service_and_data_deinterleaved_bits[i]);
// 	}


// 	memset(&dec, 0, sizeof(dec));
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
	viterbi_decoder dec2;
	dec2.decode(service_and_data_deinterleaved_bits, service_and_data_decoded_bits, bits_count+16, ntraceback);
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

// detect short preambles and feed packet decoder.
int slice_and_decode(int16_t *s, int sample_count, uint8_t *out_data, int *valid_data_len)
{
	FILE * f =fopen("ShortTraining.csv", "wb");
	fprintf(f, "N,P,A\n");

	int Nwindow = 48;

	float df_avg = 0;
	int avg_count = 0;

	bool preamble_found = false;
	float autocorrelation_throthold = 0.80;
	int preamble_start = 0;

	for(int n=0; n<sample_count-64; n++)
	{
		int64_t a[2] = {0};			// autocorrelation, [i,q]
		int64_t p = 0;				// amplitude normalizer
		int64_t _f[2] = {0};		// frequency offset, [i,q]


		for(int k=0; k<Nwindow; k++)
		{
			a[0] += s[(n+k)*2+1] * s[(n+k+16)*2+1] - s[(n+k)*2+0] * -s[(n+k+16)*2+0];
			a[1] += s[(n+k)*2+1] * -s[(n+k+16)*2+0] + s[(n+k)*2+0] * s[(n+k+16)*2+1];

			p += s[(n+k)*2+0] * s[(n+k)*2+0] + s[(n+k)*2+1] * s[(n+k)*2+1];

			//a += s[n+k] * s[n+k+16].conjugate();
			//p += s[n+k].sq_magnitude();
		}

		for(int k=0; k<16; k++)
		{
			//_f += s[n+k+Nwindow+16] * s[n+k+16+Nwindow+16].conjugate();	// note:latency ~= Nwindow+16

			_f[0] += s[(n+k+Nwindow+16)*2+1] * s[(n+k+Nwindow+32)*2+1] - s[(n+k+Nwindow+16)*2+0] * -s[(n+k+Nwindow+32)*2+0];
			_f[1] += s[(n+k+Nwindow+16)*2+1] * -s[(n+k+Nwindow+32)*2+0] + s[(n+k+Nwindow+16)*2+0] * s[(n+k+Nwindow+32)*2+1];
		}

		p = max(p, Nwindow*20);	// reject false positive in weak signal

		float autocorrelation = sqrt((float)a[0]*a[0]+a[1]*a[1]) / p;
		float df = atan2((double)_f[1], (double)_f[0])/16.0f;		// - _f.argument();

// 		fprintf(f, "%f,%f,%.0f\n", float(n), autocorrelation*1000, sqrt((float)s[(n+Nwindow+16)*2] * s[(n+Nwindow+16)*2] + s[(n+Nwindow+16)*2+1] + s[(n+Nwindow+16)*2+1]));

		if (autocorrelation > autocorrelation_throthold && !preamble_found)
		{
			df_avg += df;
			avg_count ++;

			int q = _f[1] / 65536;

			if (preamble_start == 0)
				preamble_start = n+Nwindow+16;		// note:latency ~= Nwindow+16
		}
		else if (avg_count > 3 && autocorrelation < autocorrelation_throthold)
		{
			preamble_found = true;
		}
		else
		{
			preamble_found = false;
		}

		if (avg_count > 63)
		{
			preamble_found = true;
			break;
		}
	}

	df_avg /= avg_count;

	if (preamble_found)
	{
		printf("found short training sequence @ %d sample(%fus), avg_count=%d.\n"
			"coarse frequency offset:%f degree / sample (%f Mhz)\n", preamble_start, preamble_start/20.0f, avg_count, df_avg * 180 / PI, 20 / (2 * PI / df_avg) );
	}
	else
	{
		printf("preamble not found.\n");
		return sample_count;
	}
	fclose(f);


	// copy to complex and apply coarse frequency offset;
	complex pkt[MAX_SYMBOLS_PER_SAMPLE*80+320];
	int pkt_size = min(MAX_SYMBOLS_PER_SAMPLE*80, sample_count-preamble_start);
	for(int n=0; n<pkt_size; n++)
	{
		// s2[n] = s[n] * e^(i*n*df)
		//complex offset(cos(n*df_avg), sin(n*df_avg));
		pkt[n] = complex(s[(n+preamble_start)*2+1], s[(n+preamble_start)*2+0]);	// note: reverse IQ
		pkt[n] *= cordic((n)*df_avg);//offset;
	}

// 	FILE * comp = fopen("comp_8bit.pcm", "wb");
// 	for(int i=0; i<pkt_size; i++)
// 	{
// 		int8_t I = pkt[i].real/256;
// 		int8_t Q = pkt[i].image/256;
// 		fwrite(&Q, 1, 1, comp);
// 		fwrite(&I, 1, 1, comp);
// 	}
// 	fclose(comp);


	int pkt_result = frame_decoding(pkt, pkt_size, out_data, valid_data_len);

	if (pkt_result < 0)
		return 320;

	return pkt_result + preamble_start;
}


float randf()
{
	int32_t v = ((rand()&0xff) << 24) | ((rand()&0xff) << 16) | ((rand()&0xff) << 8) | ((rand()&0xff) << 0);
	return v/2147483648.0f;
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
		int m = 512;
		int16_t I = int((s[i].real + (rand()&0xff-127)*m/float(m*2)) / m) * m;
		int16_t Q = int((s[i].image + (rand()&0xff-127)*m/float(m*2)) / m) * m;
		fwrite(&Q, 1, 2, f);
		fwrite(&I, 1, 2, f);

	}
	fclose(f);

// 	uint8_t decoded_data[4096] = {0};
// 	int decoded_count = 0;
// 	rx(s, sample_count, decoded_data, &decoded_count);

	return sample_count;
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

	fprintf(stderr, "reading....");
	char file[] = "testcases\\wifi_12mbps_1frame_8bit.pcm";
	bool _8bit = strstr(file, "8bit");
	FILE * f = fopen(file, "rb");
	if (!f)
	{
		printf("failed opening file\n");
		return 0;
	}
	fseek(f, 0, SEEK_END);
	int file_size = ftell(f);
	fseek(f, 0, SEEK_SET);

	int16_t * data;
	int sample_count;
	
	if (!_8bit)
	{
		sample_count = file_size/4;
		data = new short[sample_count*2];
		fread(data, 1, file_size, f);
	}
	else
	{
		int8_t * data8 = new int8_t[file_size];
		fread(data8, 1, file_size, f);
		sample_count = file_size/2;
		data = new int16_t[sample_count*2];
		for(int i=0; i<sample_count*2; i++)
			data[i] = data8[i];
		delete [] data8;
	}
	fclose(f);

	fprintf(stderr, "done %d samples\n", sample_count);
	int l = gettime();

	f = fopen("data.pkt", "wb");

	int pos = 0;
	int kk = 0;
	uint8_t data_out[8192];
	int valid_size = 0;
	int valid_count = 0;
	while(sample_count > pos)
	{
		int c = slice_and_decode(data+pos*2, sample_count-pos, data_out, &valid_size);
		pos += c;
		if (valid_size)
			valid_count ++ ;

		fprintf(stderr, "\r%d/%d ", valid_count, kk++);

		if (valid_size)
		{
			fwrite(&valid_size, 1, 4, f);
			fwrite(data_out, 1, valid_size, f);
		}
	}

	fclose(f);

	l = gettime() - l;

	printf("total:%d\n", l);

	return 0;
}