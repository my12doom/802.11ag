#pragma once

#include <stdint.h>
#include "complex.h"
#include "common.h"

typedef int (*mapper)(complex *symbols, uint8_t *bits);
mapper modulation2mapper(modulation modu);		// mapper: map bits to OFDM symbol
mapper modulation2demapper(modulation modu);	// demapper: map OFDM symbol to bits

int BPSK_symbols_to_bits(complex *symbols, uint8_t *bits);
int QPSK_symbols_to_bits(complex *symbols, uint8_t *bits);
int QAM16_symbols_to_bits(complex *symbols, uint8_t *bits);
int QAM64_symbols_to_bits(complex *symbols, uint8_t *bits);

int bits_to_BPSK_symbols(complex *symbols, uint8_t *bits);
int bits_to_QPSK_symbols(complex *symbols, uint8_t *bits);
int bits_to_QAM16_symbols(complex *symbols, uint8_t *bits);
int bits_to_QAM64_symbols(complex *symbols, uint8_t *bits);


