#include "stdafx.h"
#include <cmath> 
#define EPS 1.0e-10
#define JMAX 20

double qtrap2(double (*func)(double,double),double a,double b,double g)
/*改进梯形法，对于func(x,g)，当参数g已知时，在固定积分限a<x<b之间求函数func(x,g)的积分*/
{
	double trapzd2(double (*func)(double,double),double a,double b,double g,int n);
	double trapzd2();
	extern void nrerror(char error_text[]);
	int j;
	double s,olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd2(func,a,b,g,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		olds=s;
	}
	nrerror("Too many steps in routine qtrap");

	return 0.0;
}
#undef EPS
#undef JMAX

double trapzd2(double (*func)(double,double),double a,double b,double g,int n)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*((*func)(a,g)+(*func)(b,g)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += (*func)(x,g);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

