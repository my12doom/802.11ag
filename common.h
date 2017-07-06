#pragma once

#include <assert.h>

enum modulation		// values are also bits per subcarrier
{
	BPSK = 1,
	QPSK = 2,
	QAM16 = 4,
	QAM64 = 6,
};

enum puncturing
{
	_1_2 = 6,		// 6 = 1/2 * 12
	_2_3 = 8,		// 8 = 1/2 * 12
	_3_4 = 9,		// 9 = 1/2 * 12
};

static modulation rate2modulation(int mbps)
{
	if (mbps == 6)
		return BPSK;
	if (mbps == 9)
		return BPSK;
	if (mbps == 12)
		return QPSK;
	if (mbps == 18)
		return QPSK;
	if (mbps == 24)
		return QAM16;
	if (mbps == 36)
		return QAM16;
	if (mbps == 48)
		return QAM64;
	if (mbps == 54)
		return QAM64;

	assert(0);
	return BPSK;
}

static const char* modulation_name(modulation mod)
{
	if (mod == BPSK)
		return "BPSK";
	if (mod == QPSK)
		return "QPSK";
	if (mod == QAM16)
		return "16QAM";
	if (mod == QAM64)
		return "64QAM";
	return NULL;
}

static const char* puncturing_name(puncturing pun)
{
	if (pun == _1_2)
		return "1/2";
	if (pun == _2_3)
		return "2/3";
	if (pun == _3_4)
		return "3/4";
	return NULL;
}

static puncturing rate2puncturing(int mbps)
{
	if (mbps == 6)
		return _1_2;
	if (mbps == 9)
		return _3_4;
	if (mbps == 12)
		return _1_2;
	if (mbps == 18)
		return _3_4;
	if (mbps == 24)
		return _1_2;
	if (mbps == 36)
		return _3_4;
	if (mbps == 48)
		return _2_3;
	if (mbps == 54)
		return _3_4;

	assert(0);
	return _1_2;
}


static int rate_code_to_mbps(uint8_t code)
{
	if (code == 13)
		return 6;
	if (code == 15)
		return 9;
	if (code == 5)
		return 12;
	if (code == 7)
		return 18;
	if (code == 9)
		return 24;
	if (code == 11)
		return 36;
	if (code == 1)
		return 48;
	if (code == 3)
		return 54;

	return 0;
}
