#pragma once

// input: n data points
// output: a line: y=kx+b
void line_fitting(const float*x, const float *y, int point_count, float *k, float *b);

// x = 0 ~~~ point_count-1
void line_fitting_y(const float *y, int point_count, float *k, float *b);