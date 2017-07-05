#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <Windows.h>
#include <assert.h>
#include "Fourier.h"
#include "complex.h"
#include "line_fitting.h"
#include "viterbi_decoder.h"
#include "crc32_80211.h"

using namespace gr::ieee802_11;

#define countof(x) (sizeof(x)/sizeof(x[0]))
HRESULT save_bitmap(RGBQUAD *data, const wchar_t *filename, int width, int height);


int interleave_pattern[4][288];	// [BPSK, QPSK, 16QAM, 64QAM] [max: 6*48bits]

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

int rate_to_mbps(uint8_t rate)
{
	if (rate == 13)
		return 6;
	if (rate == 15)
		return 9;
	if (rate == 5)
		return 12;
	if (rate == 7)
		return 18;
	if (rate == 9)
		return 24;
	if (rate == 11)
		return 36;
	if (rate == 1)
		return 48;
	if (rate == 3)
		return 54;

	return 0;
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

complex LT_time_space[128];			// LT = long training, 2 repetition.
complex LT_frequency_space[64];		// one OFDM symbol only

int init_LTS()
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




int main()
{
	int l = GetTickCount();

	init_LTS();
	init_scrambler();
	init_interleaver_pattern();
	FILE * f = fopen("testcases\\wifi_6mbps_1frame.pcm", "rb");
	fseek(f, 0, SEEK_END);
	int file_size = ftell(f);
	fseek(f, 0, SEEK_SET);

	int sample_count = file_size/4;
	short * data = new short[file_size/2];
	fread(data, 1, file_size, f);
	fclose(f);

	complex *s = new complex[sample_count];
	for(int i=0; i<sample_count; i++)
	{
		float n = (rand()&0xff-128)*13.0f/255;
		n=0;
		s[i].real = data[i*2+1]/256 + n;		// revert IQ
		s[i].image = data[i*2+0]/256 + n;
	}

	f =fopen("ShortTraining.csv", "wb");
	fprintf(f, "N,P,A\n");

	int Nwindow = 48;

	float df_avg = 0;
	int avg_count = 0;

	bool preamble_found = false;
	float autocorrelation_throthold = 0.80;
	int preamble_start = 0;

	for(int n=0; n<sample_count-64; n++)
	{
		complex a;
		complex p;
		complex _f;
		for(int k=0; k<Nwindow; k++)
			a = a + s[n+k] * s[n+k+16].conjugate();
		for(int k=0; k<Nwindow; k++)
			p = p + s[n+k] * s[n+k].conjugate();
		for(int k=0; k<16; k++)
			_f = _f + s[n+k+Nwindow+16] * s[n+k+16+Nwindow+16].conjugate();	// note:latency ~= Nwindow+16

		float autocorrelation = a.magnitude() / p.magnitude();
		float df = _f.argument()/16.0f;
		if (autocorrelation>1)
			autocorrelation = 1;

		fprintf(f, "%f,%f,%f\n", float(n), autocorrelation*1000, s[n+Nwindow+16].magnitude());

		if (autocorrelation > autocorrelation_throthold && !preamble_found)
		{
			df_avg += df;
			avg_count ++;

// 			printf("%d:%f\n", n, df);

			if (preamble_start == 0)
				preamble_start = n+Nwindow+16;		// note:latency ~= Nwindow+16
		}
		else if (avg_count > 3 && autocorrelation < autocorrelation_throthold)
		{
			preamble_found = true;
		}
		else
		{
			printf("");
		}

		if (avg_count > 70)
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
		return -1;
	}
	fclose(f);

	// apply coarse frequency offset, write to debug file ,and swap pointer
	for(int n=0; n<sample_count; n++)
	{	
		// s2[n] = s[n] * e^(i*n*df)
		complex offset(cos(n*df_avg), sin(n*df_avg));
		s[n] = s[n] * offset;
	}

	// OFDM symbol alignment
	// use long training sequence to find first symbol
	f =fopen("alignment.csv", "wb");
	fprintf(f, "N,P,A\n");
	float autocorrection_value[3] = {0};
	int autocorrection_pos[3] = {0};
	for(int n=80; n<320+50; n++)			// 160: number of sample of short training sequence
		// 320: number of sample of long training sequence
	{
		complex a;
		float p = 0;
		for(int k=0; k<64; k++)
		{
			a = a + s[preamble_start+n+k] * LT_time_space[k].conjugate();
			p += (s[preamble_start+n+k] * s[preamble_start+n+k].conjugate()).real;
		}

		float autocorelation = a.magnitude()/p;
		fprintf(f, "%d,%f\n", n, autocorelation);

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

// 	peak_pos += 64;
	printf("long training sequence @ %d of preamble\n", peak_pos);
	int symbol_start = preamble_start + peak_pos + 64;
	printf("symbols start @ %d (%.1fus)\n", symbol_start, symbol_start/20.0f);

	// fine frequency offset
	complex _df;
	float sigma = 0;
	for(int i=0; i<64; i++)
	{
		float v_df = 0.005*PI/180;
		complex offset1(cos((i)*v_df), sin((i)*v_df));
		complex offset2(cos((i+64)*v_df), sin((i+64)*v_df));
		complex _f = s[symbol_start-64+i] * s[symbol_start-128+i].conjugate();
		_df = _df + _f;

		float _f2 = phase_sub(s[symbol_start-64+i].argument(), s[symbol_start-128+i].argument());
// 		printf("%f\t%f\n", _f2*180/PI, _f.argument()*180/PI);

		sigma += _f2;
	}

	float df = _df.argument()/64;
// 	sigma /= 64;
// 
// 	df = sigma * PI / 180;

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

		assert(abs(phase_diffs- (LT_frequency_space[n]*lt_rx_fft[n].conjugate()).argument()) < 0.05);
		complex normalized = lt_rx_fft[n]*h[n];

// 		printf("%d,%.2f,%.4f\n", i, normalized.argument(), normalized.magnitude());
	}


	// "SIGNAL" symbol, BPSK, 1/2 convolutional coded
	// if this symbol fails to decode, we drop the whole frame
	complex signal_fft[64];
	fft_complex(&s[symbol_start+16], signal_fft);
	uint8_t signal_bits[48] = {0};
	uint8_t signal_bits_deinterleaved[64] = {0};
	uint8_t signal_decoded_bits[48];
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

	int rate = (signal_decoded_bits[0] << 3) | (signal_decoded_bits[1] << 2) | (signal_decoded_bits[2] << 1) | signal_decoded_bits[3];

	printf("rate:%dMbps\n", rate_to_mbps(rate));
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

	if (rate_to_mbps(rate) == 0
		|| signal_decoded_bits[4] != 0
		|| parity != signal_decoded_bits[17]
		|| !tail_is_zero
		)
	{
		printf("invalid signal field\n");
		return -1;
	}


	int data_symbol_count = (length*8+ 16+6 + rate_to_mbps(rate)*4-1) / (rate_to_mbps(rate)*4);
	printf("%d data symbols\n", data_symbol_count);

	int symbol_count = data_symbol_count + 1;		// 1 = signal

	// fft all symbols
	complex symbols[200][64];
	for(int i=0; i<(sample_count-symbol_start)/80; i++)
	{
		fft_complex(s+symbol_start+i*80+16, symbols[i]);

		// equalizing
		float phase_offset_per_symbol = -df * 80;
		float phase_offset = (i+1)*phase_offset_per_symbol;
		complex offset(cos(phase_offset), sin(phase_offset));

		for(int j=0; j<64; j++)
			symbols[i][j] = symbols[i][j] * h[j] * offset;
	}


	// compensate for sampling rate error, which introduce a triangular shaped phase shift to all subcarriers
	// use pilot tones to compensate for it
	float pilot_phase[4][200];
	float pilot_phase_error[4][200];
	for(int i=0; i<(sample_count-symbol_start)/80; i++)
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
			phase_error[j][idx] = k_i * j;
	}

	FILE * csv = fopen("constellation.csv", "wb");
	fprintf(f, "N,P,A\n");

	uint8_t service_and_data_bits[48*512];
	uint8_t service_and_data_deinterleaved_bits[48*512];
	uint8_t service_and_data_decoded_bits[24*512];
	int pos = 0;
	for(int i=1; i<symbol_count; i++)
	{
		for(int j=-26; j<=26; j++)
		{
			if (j == 0)
				continue;

			int n = j>0?j:j+64;
			if ((n==7||n==21||n==-7||n==-21||n==64-7||n==64-21))
				continue;

			complex offset2(cos(phase_error[i][n]), -sin(phase_error[i][n]));
			complex v = symbols[i][n] * offset2;

			fprintf(csv, "%d,%f,%f,%f\n", j, v.real, v.image, v.argument());
			service_and_data_bits[pos++] = map_BPSK(v);
		}

		printf("");
	}

	fclose(csv);

	for(int i=0; i<48*data_symbol_count; i++)
	{
		int symbol = i/48;
		service_and_data_deinterleaved_bits[i] = service_and_data_bits[i/48*48 + interleave_pattern[0][i%48]];
	}

	dec.decode(service_and_data_deinterleaved_bits, service_and_data_decoded_bits, 48*data_symbol_count);
	uint8_t out_bytes[4096];
	descramble(service_and_data_decoded_bits, out_bytes, 24*data_symbol_count);

	uint32_t crc = crc32_80211(out_bytes+2, length-4);

	printf("FCS(calculated)=0x%08x\n", crc);
	printf("FCS(received)=0x%02x%02x%02x%02x\n", out_bytes[2+length-4+3], out_bytes[2+length-4+2], out_bytes[2+length-4+1], out_bytes[2+length-4+0]);

	f = fopen("data.bin", "wb");
	fwrite(out_bytes+2, 1, length, f);
	fclose(f);

	printf("%dms\n", GetTickCount()-l);
	return 0;
}
