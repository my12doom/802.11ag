#include "line_fitting.h"

void line_fitting(const float*x, const float *y, int point_count, float *k, float *b)
{
	float xsq = 0;
	float xy = 0;
	float _x = 0;
	float _y = 0;

	for(int i=0; i<point_count; i++)
	{
		xsq += x[i]*x[i];
		xy += x[i]*y[i];
		_x += x[i];
		_y += y[i];
	}

	*k = (point_count *xy - _x * _y) / (point_count*xsq - _x*_x);
	*b = (xsq*_y - _x*xy) / (point_count*xsq - _x*_x);
}

void line_fitting_y(const float *y, int point_count, float *k, float *b)
{
	static float x[1024] = {-1};
	if (x[0] == -1)
	{
		for(int i=0; i<sizeof(x)/sizeof(float); i++)
			x[i] = i;
	}

	line_fitting(x, y, point_count, k, b);
}