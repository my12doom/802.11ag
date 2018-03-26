#include "convolutional.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>

#define MAX_BITS (4096*8)



static uint8_t ones[256];
static int16_t polyA[128];
static int16_t polyB[128];
static int16_t scale;
static int metric[MAX_BITS][64];
static uint16_t path[MAX_BITS];

static int imin(int a, int b)
{
	return a>b?b:a;
}

int convolutional80211_init(int16_t scale)
{
	::scale = scale;

	for(int n=0; n<256; n++)
	{
		ones[n] = 0;
		for(int i = 0; i < 8; i++) 
		{
			if(n & (1 << i)) 
			{
				ones[n]++;
			}
		}
	}
	
	for(int n=0; n<128; n++)
	{
		polyA[n] = (ones[n & 0x6d] & 1) ? scale : -scale;	// 802.11 generator polynomial A: binary 1101101(0x6d, 0155)
		polyB[n] = (ones[n & 0x4f] & 1) ? scale : -scale;	// 802.11 generator polynomial B: binary 1001111(0x4f, 0117)
	}

	return 0;
}

int convolutional80211_encode(const uint8_t *input_bits, int16_t *output, int input_count)
{
	uint8_t state = 0;
	for(int i = 0; i < input_count; i++) 
	{
		assert(input_bits[i] == 0 || input_bits[i] == 1);
		state = ((state&0x3f) << 1) | input_bits[i];			// shift left and input data on LSB
		*output++ = polyA[state];
		*output++ = polyB[state];
	}

	return 0;
}

int viterbi_soft_decode(const int16_t *inputs, int input_count, uint8_t *output)
{
	if (input_count >= MAX_BITS*2)
		return -1;
	int output_count = input_count/2;

	// build trellis, note that we ignore starting and ending part, leaving that part to trace back stage
	memset(metric, 0, sizeof(metric));

	for(int i=0; i<output_count; i++)
	{
		int16_t inA = inputs[(i<<1)];
		int16_t inB = inputs[(i<<1)|1];

		for(int j=0; j<64; j++)
		{
			int16_t state0 = j>>1;
			int16_t state1 = (j>>1) | (1<<5);
			int16_t LSB = j&1;

			int metric_left0 = i ? metric[i-1][state0] : 0;
			int metric_left1 = i ? metric[i-1][state1] : 0;
			
			int16_t out0[2] = {polyA[(state0<<1)|LSB], polyB[(state0<<1)|LSB]};
			int16_t out1[2] = {polyA[(state1<<1)|LSB], polyB[(state1<<1)|LSB]};

			int metric0 = metric_left0 + (out0[0]-inA)*(out0[0]-inA) + (out0[1]-inB)*(out0[1]-inB);
			int metric1 = metric_left1 + (out1[0]-inA)*(out1[0]-inA) + (out1[1]-inB)*(out1[1]-inB);

			metric[i][j] = imin(metric0, metric1);
		}
	}

	// trace back
	int least_metric = 0x7fffffff;
	int least_index;
	for(int j=0; j<64; j++)
	{
		if (metric[output_count-1][j] < least_metric)
		{
			least_index = j;
			least_metric = metric[output_count-1][j];
		}
	}

	for(int i=output_count-1; i>0; i--)
	{
		output[i] = least_index & 1;

		int16_t pre0 = least_index>>1;
		int16_t pre1 = (least_index>>1) | (1<<5);

		// clamp unreachable starting phase of trellis
		if (i<=6)
			pre1 = pre0;

		least_index = metric[i-1][pre0] < metric[i-1][pre1] ? pre0 : pre1;
	}

	output[0] = least_index & 1;

	return least_metric;
}

int convolutional80211_test()
{
	convolutional80211_init();

	const int count = 10;
	uint8_t bits[count];
	for(int i=0; i<count; i++)
		bits[i] = (i >> 1)&1;

	int16_t encoded[count*2];
	convolutional80211_encode(bits, encoded, count);
	for(int i=1; i<count*2; i+=2)
		encoded[i] = 0;

	uint8_t decoded[count];
	viterbi_soft_decode(encoded, count*2, decoded);

	int bit_error = 0;
	for(int i=0; i<count; i++)
		if (bits[i] != decoded[i])
			bit_error++;
	return 0;
}
