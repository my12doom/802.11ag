#pragma once

#include <stdint.h>

int convolutional80211_init(int16_t scale = 127);
int convolutional80211_encode(const uint8_t *input_bits, int16_t *output, int input_count);

// use 0 for punctured bits
int viterbi_soft_decode(const int16_t *inputs, int input_count, uint8_t *output);


int convolutional80211_test();