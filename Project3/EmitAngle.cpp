#include "stdafx.h"
#include <cmath>
#define	 PI          3.1415926535
#define  ONE_MINUS_COSZERO 1.0E-12

void EmitAngle(double Angle[], int NA, double x, double y, double z, double w)
{
	double A, B, C, D;
	double L1,L2,L3,L4;
	double d1,d2,d3,d4;
	double w2;
	double ZA_max, ZA_min, ZA, AA_max, AA_min, AA;
	double Tmax, Tmin;
	double cosmax,cosmin, sinmax,sinmin;
	int    i;

	w2 = w/2.0;
	A = atan(fabs(w2-y)/fabs(w2-x));
	B = atan(fabs(w2-y)/fabs(-w2-x));
	C = atan(fabs(-w2-y)/fabs(-w2-x));
	D = atan(fabs(-w2-y)/fabs(w2-x));

	L1 = sqrt(pow(x-w2,2) + pow(y-w2,2));
	L2 = sqrt(pow(x+w2,2) + pow(y-w2,2));
	L3 = sqrt(pow(x+w2,2) + pow(y+w2,2));
	L4 = sqrt(pow(x-w2,2) + pow(y+w2,2));

	d1 = fabs(x-w2);
	d2 = fabs(x+w2);
	d3 = fabs(y-w2);
	d4 = fabs(y+w2);

	Tmax = L1;
	if ( L2>Tmax ) Tmax = L2;
	if ( L3>Tmax ) Tmax = L3;
	if ( L4>Tmax ) Tmax = L4;

	if ( fabs(x)>w2 && fabs(y)>w2 ) {
		Tmin = L1;
		if ( L2<Tmin ) Tmin = L2;
		if ( L3<Tmin ) Tmin = L3;
		if ( L4<Tmin ) Tmin = L4;
	} else if ( fabs(x)>w2 || fabs(y)>w2 ) {
		if ( fabs(x)>w2 ) {
			Tmin = d1;
			if (d2<Tmin) Tmin = d2;
		} else {
			Tmin = d3;
			if (d4<Tmin) Tmin = d4;
		}
	} else {
		Tmin = 0.0;
		Tmax = d1;
		if ( d2<Tmax ) Tmax= d2;
		if ( d3<Tmax ) Tmax = d3;
		if ( d4<Tmax ) Tmax = d4;
	}
//	cout << "Tmax= " << Tmax << "    Tmin= " << Tmin << endl;
	ZA_max = atan(Tmax/z);
	ZA_min = atan(Tmin/z);

	if ( fabs(x-w2)<ONE_MINUS_COSZERO ) {
		if ( y>(-w2) && y<w2) {
			AA_min = 0.5 * PI;
			AA_max = 1.5 * PI;
		} else if ( y>=w2 ) {
			AA_min = PI + B;
			AA_max = 1.5 * PI;
		} else {
			AA_min = 0.5 * PI;
			AA_max = PI - C;
		}
	} else if ( fabs(x+w2)<ONE_MINUS_COSZERO ) {
		if ( y>(-w2) && y<w2) {
			AA_min = -0.5 * PI;
			AA_max = 0.5 * PI;
		} else if ( y>=w2 ) {
			AA_min = 1.5 * PI;
			AA_max = 2 * PI - A;
		} else {
			AA_min = D;
			AA_max = 0.5 * PI;
		}
	} else if ( x>w2 ) {
		if ( y>(-w2) && y<w2) {
			AA_min = PI - A;
			AA_max = PI + D;
		} else if ( y>=w2 ) {
			AA_min = PI + B;
			AA_max = PI + D;
		} else {
			AA_min = PI - A;
			AA_max = PI - C;
		}
	} else if ( x<(-w2)) {
		if ( y>(-w2) && y<w2) {
			AA_min = -C;
			AA_max = B;
		} else if ( y>=w2 ) {
			AA_min = 2 * PI - C;
			AA_max = 2 * PI - A;
		} else {
			AA_min = D;
			AA_max = B;
		}
	} else {
		if ( y>(-w2) && y<w2) {
			AA_min = 0;
			AA_max = 2 * PI;
		} else if ( y>=w2 ) {
			AA_min = PI + B;
			AA_max = 2 * PI - A;
		} else {
			AA_min = D;
			AA_max = PI - C;
		}
	}
	Angle[0] = ZA_max;
	Angle[1] = ZA_min;
	Angle[2] = AA_max;
	Angle[3] = AA_min;

}
#undef	 PI 
#undef   ONE_MINUS_COSZERO
