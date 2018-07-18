#ifndef FIX_FFT_H
#define FIX_FFT_H

int fix_fft(short fr[], short fi[], short m, short inverse);
void fft_fixed(unsigned int p_nSamples, bool p_bInverseTransform, float *p_lpRealIn, float *p_lpImagIn, float *p_lpRealOut, float *p_lpImagOut);

#endif