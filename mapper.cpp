#include "mapper.h"

mapper modulation2mapper(modulation modu)		// mapper: map bits to OFDM symbol
{
	if (modu == BPSK)
		return bits_to_BPSK_symbols;
	if (modu == QPSK)
		return bits_to_QPSK_symbols;
	if (modu == QAM16)
		return bits_to_QAM16_symbols;
	if (modu == QAM64)
		return bits_to_QAM64_symbols;

	return NULL;
}
mapper modulation2demapper(modulation modu)	// demapper: map OFDM symbol to bits
{
	if (modu == BPSK)
		return BPSK_symbols_to_bits;
	if (modu == QPSK)
		return QPSK_symbols_to_bits;
	if (modu == QAM16)
		return QAM16_symbols_to_bits;
	if (modu == QAM64)
		return QAM64_symbols_to_bits;

	return NULL;
}

int BPSK_symbols_to_bits(complex *symbols, uint8_t *bits)
{
	int j = 0;
	for(int i=-26;i<=26;i++)
	{
		if (i==0 || i == -7 || i == -21 || i == 7 || i == 21)
			continue;
		int idx = i>0?i:i+64;
		bits[j++] = symbols[idx].real > 0 ? 1 : 0;
	}

	assert(j == 48);

	return j;
}

int QPSK_symbols_to_bits(complex *symbols, uint8_t *bits)
{
	int j = 0;
	for(int i=-26;i<=26;i++)
	{
		if (i==0 || i == -7 || i == -21 || i == 7 || i == 21)
			continue;
		int idx = i>0?i:i+64;
		int b0 = symbols[idx].real > 0 ? 1 : 0;
		int b1 = symbols[idx].image < 0 ? 1 : 0;	// WTF? image reversed?

		bits[j++] = b0;
		bits[j++] = b1;

	}

	assert(j == 96);

	return j;
}

int QAM16_symbols_to_bits(complex *symbols, uint8_t *bits)
{
	int j = 0;
	for(int i=-26;i<=26;i++)
	{
		if (i==0 || i == -7 || i == -21 || i == 7 || i == 21)
			continue;
		int idx = i>0?i:i+64;

		int b0 = symbols[idx].real > 0 ? 1 : 0;
		int b1 = abs(symbols[idx].real) < 0.666f ? 1 : 0;
		int b2 = symbols[idx].image < 0 ? 1 : 0;			// WTF? image reversed?
		int b3 = abs(symbols[idx].image) < 0.666f ? 1 : 0;

		bits[j++] = b0;
		bits[j++] = b1;
		bits[j++] = b2;
		bits[j++] = b3;
	}

	assert(j == 192);
	return j;
}

int QAM64_symbols_to_bits(complex *symbols, uint8_t *bits)
{
	int j = 0;
	for(int i=-26;i<=26;i++)
	{
		if (i==0 || i == -7 || i == -21 || i == 7 || i == 21)
			continue;
		int idx = i>0?i:i+64;

		int b0 = symbols[idx].real > 0 ? 1 : 0;
		int b1 = abs(symbols[idx].real) < 0.5715f ? 1 : 0;
		int b2 = (abs(symbols[idx].real) > 0.2857f && abs(symbols[idx].real) < 0.8571) ? 1 : 0;

		int b3 = symbols[idx].image < 0 ? 1 : 0;	// WTF? image reversed?
		int b4 = abs(symbols[idx].image) < 0.5715f ? 1 : 0;
		int b5 = (abs(symbols[idx].image) > 0.2857f && abs(symbols[idx].image) < 0.8571) ? 1 : 0;

		bits[j++] = b0;
		bits[j++] = b1;
		bits[j++] = b2;
		bits[j++] = b3;
		bits[j++] = b4;
		bits[j++] = b5;
	}

	assert(j == 288);
	return j;
}

int bits_to_BPSK_symbols(complex *symbols, uint8_t *bits)
{
	return 0;
}
int bits_to_QPSK_symbols(complex *symbols, uint8_t *bits)
{
	return 0;
}
int bits_to_QAM16_symbols(complex *symbols, uint8_t *bits)
{
	return 0;
}
int bits_to_QAM64_symbols(complex *symbols, uint8_t *bits)
{
	return 0;
}
