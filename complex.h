#pragma once

class complex
{
public:
	complex()
	{
		this->real = 0;
		this->image = 0;
	}
	complex(float real, float image)
	{
		this->real = real;
		this->image = image;
	}
	complex(const complex &v)
	{
		this->real = v.real;
		this->image = v.image;
	}
	static complex from_phase_magnitude(float phase, float magnitude)
	{
		complex o;
		o.real = cos(phase) * magnitude;
		o.image = sin(phase) * magnitude;
		return o;
	}
	~complex()
	{

	}

	float real;
	float image;
	float sq_magnitude() 
	{
		return real*real+image*image;
	}
	float magnitude() 
	{
		return sqrt(real*real+image*image);
	}
	float argument() 
	{
		return atan2(image, real);
	}
	complex conjugate() 
	{
		complex o(real, -image);
		return o;
	}

	complex operator *(const complex &b)
	{
		complex o;

		o.real = this->real * b.real - this->image * b.image;
		o.image = this->real * b.image + this->image * b.real;

		return o;
	}

	complex operator *=(const complex &b)
	{
		complex o;

		o.real = this->real * b.real - this->image * b.image;
		o.image = this->real * b.image + this->image * b.real;

		*this = o;
		return o;
	}

	complex operator +(const complex &b)
	{
		complex o(real, image);

		o.real += b.real;
		o.image += b.image;

		return o;
	}

	complex operator /(const float &b)
	{
		complex o(real, image);

		o.real /= b;
		o.image /= b;

		return o;
	}

	complex operator =(const float &b)
	{
		this->real = b;
		this->image = 0;

		return *this;
	}
};
