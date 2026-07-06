#include "stdafx.h"
#include <cmath> 
#define EPS 1.0e-10
#define JMAX 20

double qtrap(double (*func)(double),double a,double b)
/*改进梯形法，在固定积分限a和b之间求函数func(x)的积分*/
{
	double trapzd(double (*func)(double),double a,double b,int n);
	double trapzd();
	extern void nrerror(char error_text[]);
	int j;
	double s,olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd(func,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		olds=s;
	}
	nrerror("Too many steps in routine qtrap");

	return 0.0;
}
#undef EPS
#undef JMAX

double trapzd(double (*func)(double),double a,double b,int n)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*((*func)(a)+(*func)(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += (*func)(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

