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
		int b1 = fabs(symbols[idx].real) < 0.666f ? 1 : 0;
		int b2 = symbols[idx].image < 0 ? 1 : 0;			// WTF? image reversed?
		int b3 = fabs(symbols[idx].image) < 0.666f ? 1 : 0;

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
		int b1 = fabs(symbols[idx].real) < 0.5715f ? 1 : 0;
		int b2 = (fabs(symbols[idx].real) > 0.2857f && fabs(symbols[idx].real) < 0.8571) ? 1 : 0;

		int b3 = symbols[idx].image < 0 ? 1 : 0;	// WTF? image reversed?
		int b4 = fabs(symbols[idx].image) < 0.5715f ? 1 : 0;
		int b5 = (fabs(symbols[idx].image) > 0.2857f && fabs(symbols[idx].image) < 0.8571) ? 1 : 0;

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
	static complex tbl[4] =
	{
		complex(-1.0f, 1.0f),			// WTF? image reversed?
		complex(-1.0f, -1.0f),
		complex(1.0f, 1.0f),
		complex(1.0f, -1.0f),
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
	static complex tbl[16] =
	{
		complex(-1.0f, 1.0f),			// WTF? image reversed?
		complex(-1.0f, 0.3333f),
		complex(-1.0f, -1.0f),
		complex(-1.0f, -0.3333f),

		complex(-0.3333f, 1.0f),			// WTF? image reversed?
		complex(-0.3333f, 0.3333f),
		complex(-0.3333f, -1.0f),
		complex(-0.3333f, -0.3333f),

		complex(1.0f, 1.0f),			// WTF? image reversed?
		complex(1.0f, 0.3333f),
		complex(1.0f, -1.0f),
		complex(1.0f, -0.3333f),

		complex(0.3333f, 1.0f),			// WTF? image reversed?
		complex(0.3333f, 0.3333f),
		complex(0.3333f, -1.0f),
		complex(0.3333f, -0.3333f),
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
	int j = 0;
	static complex tbl[64] =
	{
		complex(-1.0f, 1.0f),
		complex(-1.0f, 0.7143f),
		complex(-1.0f, 0.1429f),
		complex(-1.0f, 0.4286f),
		complex(-1.0f, -1.0f),
		complex(-1.0f, -0.7143f),
		complex(-1.0f, -0.1429f),
		complex(-1.0f, -0.4286f),

		complex(-0.7143f, 1.0f),
		complex(-0.7143f, 0.7143f),
		complex(-0.7143f, 0.1429f),
		complex(-0.7143f, 0.4286f),
		complex(-0.7143f, -1.0f),
		complex(-0.7143f, -0.7143f),
		complex(-0.7143f, -0.1429f),
		complex(-0.7143f, -0.4286f),

		complex(-0.1429f, 1.0f),
		complex(-0.1429f, 0.7143f),
		complex(-0.1429f, 0.1429f),
		complex(-0.1429f, 0.4286f),
		complex(-0.1429f, -1.0f),
		complex(-0.1429f, -0.7143f),
		complex(-0.1429f, -0.1429f),
		complex(-0.1429f, -0.4286f),

		complex(-0.4286f, 1.0f),
		complex(-0.4286f, 0.7143f),
		complex(-0.4286f, 0.1429f),
		complex(-0.4286f, 0.4286f),
		complex(-0.4286f, -1.0f),
		complex(-0.4286f, -0.7143f),
		complex(-0.4286f, -0.1429f),
		complex(-0.4286f, -0.4286f),

		//
		complex(+1.0f, 1.0f),
		complex(+1.0f, 0.7143f),
		complex(+1.0f, 0.1429f),
		complex(+1.0f, 0.4286f),
		complex(+1.0f, -1.0f),
		complex(+1.0f, -0.7143f),
		complex(+1.0f, -0.1429f),
		complex(+1.0f, -0.4286f),

		complex(+0.7143f, 1.0f),
		complex(+0.7143f, 0.7143f),
		complex(+0.7143f, 0.1429f),
		complex(+0.7143f, 0.4286f),
		complex(+0.7143f, -1.0f),
		complex(+0.7143f, -0.7143f),
		complex(+0.7143f, -0.1429f),
		complex(+0.7143f, -0.4286f),

		complex(+0.1429f, 1.0f),
		complex(+0.1429f, 0.7143f),
		complex(+0.1429f, 0.1429f),
		complex(+0.1429f, 0.4286f),
		complex(+0.1429f, -1.0f),
		complex(+0.1429f, -0.7143f),
		complex(+0.1429f, -0.1429f),
		complex(+0.1429f, -0.4286f),

		complex(+0.4286f, 1.0f),
		complex(+0.4286f, 0.7143f),
		complex(+0.4286f, 0.1429f),
		complex(+0.4286f, 0.4286f),
		complex(+0.4286f, -1.0f),
		complex(+0.4286f, -0.7143f),
		complex(+0.4286f, -0.1429f),
		complex(+0.4286f, -0.4286f),

	};

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
