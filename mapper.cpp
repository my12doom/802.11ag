#include "mapper.h"

mapper get_mapper(modulation modu)		// mapper: map bits to OFDM symbol
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
mapper get_demapper(modulation modu)	// demapper: map OFDM symbol to bits
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
	float level = 2.0 * sqrt(0.1);
	for(int i=-26;i<=26;i++)
	{
		if (i==0 || i == -7 || i == -21 || i == 7 || i == 21)
			continue;
		int idx = i>0?i:i+64;

		int b0 = symbols[idx].real > 0 ? 1 : 0;
		int b1 = fabs(symbols[idx].real) < level ? 1 : 0;
		int b2 = symbols[idx].image < 0 ? 1 : 0;			// WTF? image reversed?
		int b3 = fabs(symbols[idx].image) < level ? 1 : 0;

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
	float level1 = 2/sqrt(42.0);
	float level2 = 4/sqrt(42.0);
	float level3 = 6/sqrt(42.0);
	int j = 0;
	for(int i=-26;i<=26;i++)
	{
		if (i==0 || i == -7 || i == -21 || i == 7 || i == 21)
			continue;
		int idx = i>0?i:i+64;

		int b0 = symbols[idx].real > 0 ? 1 : 0;
		int b1 = fabs(symbols[idx].real) < level2 ? 1 : 0;
		int b2 = (fabs(symbols[idx].real) > level1 && fabs(symbols[idx].real) < level3) ? 1 : 0;

		int b3 = symbols[idx].image < 0 ? 1 : 0;	// WTF? image reversed?
		int b4 = fabs(symbols[idx].image) < level2 ? 1 : 0;
		int b5 = (fabs(symbols[idx].image) > level1 && fabs(symbols[idx].image) < level3) ? 1 : 0;

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
	int j = 0;
	for(int i=-26;i<=26;i++)
	{
		if (i==0 || i == -7 || i == -21 || i == 7 || i == 21)
			continue;
		int idx = i>0?i:i+64;
		symbols[idx].real = bits[j++] ? 1.0f : -1.0f;
	}

	assert(j == 48);

	return j;
}

int bits_to_QPSK_symbols(complex *symbols, uint8_t *bits)
{
	int j = 0;
	static float level = sqrt(2.0)/2;
	static complex tbl[4] =
	{
		complex(-level, level),			// WTF? image reversed?
		complex(-level, -level),
		complex(level, level),
		complex(level, -level),
	};
	for(int i=-26;i<=26;i++)
	{
		if (i==0 || i == -7 || i == -21 || i == 7 || i == 21)
			continue;
		int idx = i>0?i:i+64;
		int symbol = (bits[j] << 1) | bits[j+1];
		j += 2;
		symbols[idx] = tbl[symbol];
	}

	assert(j == 96);

	return j;
}

int bits_to_QAM16_symbols(complex *symbols, uint8_t *bits)
{
	int j = 0;
	static float level = sqrt(0.1);
	static complex tbl[16] =
	{
		complex(-3*level, 3*level),			// WTF? image reversed?
		complex(-3*level, 1*level),
		complex(-3*level, -3*level),
		complex(-3*level, -1*level),

		complex(-1*level, 3*level),			// WTF? image reversed?
		complex(-1*level, 1*level),
		complex(-1*level, -3*level),
		complex(-1*level, -1*level),

		complex(3*level, 3*level),			// WTF? image reversed?
		complex(3*level, 1*level),
		complex(3*level, -3*level),
		complex(3*level, -1*level),

		complex(1*level, 3*level),			// WTF? image reversed?
		complex(1*level, 1*level),
		complex(1*level, -3*level),
		complex(1*level, -1*level),
	};

	for(int i=-26;i<=26;i++)
	{
		if (i==0 || i == -7 || i == -21 || i == 7 || i == 21)
			continue;
		int idx = i>0?i:i+64;
		int symbol = (bits[j] << 3) | (bits[j+1] << 2) | (bits[j+2] << 1) | (bits[j+3] << 0);
		j += 4;
		symbols[idx] = tbl[symbol];
	}

	assert(j == 192);

	return j;
}

int bits_to_QAM64_symbols(complex *symbols, uint8_t *bits)
{
	static float level = sqrt(1/42.0);
	static complex tbl[64] =
	{
		complex(-7 * level, 7 * level),
		complex(-7 * level, 5 * level),
		complex(-7 * level, 1 * level),
		complex(-7 * level, 3 * level),
		complex(-7 * level, -7 * level),
		complex(-7 * level, -5 * level),
		complex(-7 * level, -1 * level),
		complex(-7 * level, -3 * level),

		complex(-5 * level, 7 * level),
		complex(-5 * level, 5 * level),
		complex(-5 * level, 1 * level),
		complex(-5 * level, 3 * level),
		complex(-5 * level, -7 * level),
		complex(-5 * level, -5 * level),
		complex(-5 * level, -1 * level),
		complex(-5 * level, -3 * level),

		complex(-1 * level, 7 * level),
		complex(-1 * level, 5 * level),
		complex(-1 * level, 1 * level),
		complex(-1 * level, 3 * level),
		complex(-1 * level, -7 * level),
		complex(-1 * level, -5 * level),
		complex(-1 * level, -1 * level),
		complex(-1 * level, -3 * level),

		complex(-3 * level, 7 * level),
		complex(-3 * level, 5 * level),
		complex(-3 * level, 1 * level),
		complex(-3 * level, 3 * level),
		complex(-3 * level, -7 * level),
		complex(-3 * level, -5 * level),
		complex(-3 * level, -1 * level),
		complex(-3 * level, -3 * level),

		//
		complex(+7 * level, 7 * level),
		complex(+7 * level, 5 * level),
		complex(+7 * level, 1 * level),
		complex(+7 * level, 3 * level),
		complex(+7 * level, -7 * level),
		complex(+7 * level, -5 * level),
		complex(+7 * level, -1 * level),
		complex(+7 * level, -3 * level),

		complex(+5 * level, 7 * level),
		complex(+5 * level, 5 * level),
		complex(+5 * level, 1 * level),
		complex(+5 * level, 3 * level),
		complex(+5 * level, -7 * level),
		complex(+5 * level, -5 * level),
		complex(+5 * level, -1 * level),
		complex(+5 * level, -3 * level),

		complex(+1 * level, 7 * level),
		complex(+1 * level, 5 * level),
		complex(+1 * level, 1 * level),
		complex(+1 * level, 3 * level),
		complex(+1 * level, -7 * level),
		complex(+1 * level, -5 * level),
		complex(+1 * level, -1 * level),
		complex(+1 * level, -3 * level),

		complex(+3 * level, 7 * level),
		complex(+3 * level, 5 * level),
		complex(+3 * level, 1 * level),
		complex(+3 * level, 3 * level),
		complex(+3 * level, -7 * level),
		complex(+3 * level, -5 * level),
		complex(+3 * level, -1 * level),
		complex(+3 * level, -3 * level),

	};

	int j = 0;
	for(int i=-26;i<=26;i++)
	{
		if (i==0 || i == -7 || i == -21 || i == 7 || i == 21)
			continue;
		int idx = i>0?i:i+64;
		int symbol = (bits[j] << 5) | (bits[j+1] << 4) | (bits[j+2] << 3) | 
				(bits[j+3] << 2)| (bits[j+4] << 1) | (bits[j+5] << 0);
		j += 6;
		symbols[idx] = tbl[symbol];
	}

	assert(j == 288);

	return j;
}
