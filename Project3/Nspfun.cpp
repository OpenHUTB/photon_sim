#include "stdafx.h"
#include <cmath>
#define  PI  3.1415926535

void Nspfun(char typespf, int N, double dataX[], double dataY[], double G)
/*归一化的散射相函数*/
{
	extern double qtrap(double (*func)(double),double a,double b);
	extern double qtrap2(double (*func)(double,double),double a,double b,double g);
	double Rspfun(double x);
	double HGspfun(double x,double g);
	double IHGspfun(double x,double g);
	ofstream SPFM, SPFA;
	double a, b, s, base, u;
	a = 0.0;

	switch(typespf){
		case 'R':{
			SPFM.open("Normalized Molecule Scatter Phase Function.dat");
	        if (!SPFM){
		        cerr << "error:unable to open Scatter Phase Function file" << endl;
		        exit(1);
	        }
	        SPFM.setf(ios_base::scientific, ios_base::floatfield);
			SPFM << "Scatter phase function type: Molecule, Rayleigh." << endl;

			base = qtrap(Rspfun,a,PI);		// 梯形的面积
			for(int i=1;i<=N;i++){
				b = i*PI/N;
				s = qtrap(Rspfun,a,b); 
				u = s/base;
				dataX[i-1] = b;				// 记录下x，一直到Pi
				dataY[i-1] = u;				// 记录下一小片面积在整个梯形面积的百分比
				SPFM.width(15);
				SPFM << dataX[i-1];			// 写入文件：Normalized Molecule Scatter Phase Function.dat
				SPFM.width(15);
				SPFM << dataY[i-1] << endl;
			}
			break;
				 }
		case 'H':{
			SPFA.open("Normalized Aerosol Scatter Phase Function.dat");
	        if (!SPFA){
		        cerr << "error:unable to open Scatter Phase Function file" << endl;
		        exit(1);
	        }      
			SPFA << "Scatter phase function type: Aerosol, H-G." << endl;
			SPFA << "asymetry factor: " << G << endl;
			SPFA.setf(ios_base::scientific, ios_base::floatfield);

			base = qtrap2(HGspfun,a,PI,G);
			for(int i=1;i<=N;i++){
				b = i*PI/N;
			    s = qtrap2(HGspfun,a,b,G); 
			    u = s/base;
				dataX[i-1] = b;
				dataY[i-1] = u;
				SPFA.width(15);
				SPFA << dataX[i-1];
				SPFA.width(15);
				SPFA << dataY[i-1] << endl;
			}
			break;
				 }
		case 'I':{
			SPFA.open("Normalized Aerosol Scatter Phase Function.dat");
	        if (!SPFA){
		        cerr << "error:unable to open Scatter Phase Function file" << endl;
		        exit(1);
	        }
			SPFA << "Scatter phase function type: Aerosol, Improved H-G." << endl;
			SPFA << "asymetry factor: " << G << endl;
	        SPFA.setf(ios_base::scientific, ios_base::floatfield);

			base = qtrap2(IHGspfun,a,PI,G);
			for(int i=1;i<=N;i++){
			    b = i*PI/N;
			    s = qtrap2(IHGspfun,a,b,G); 
			    u = s/base;
				dataX[i-1] = b;
				dataY[i-1] = u;				
				SPFA.width(15);
				SPFA << dataX[i-1];
				SPFA.width(15);
				SPFA << dataY[i-1] << endl;
			}
			break;
				 }
		default:cout<<"scatter phase function is wrong!" << endl;
	}
	SPFM.close();
	SPFA.close();
}

double Rspfun(double x)
/*瑞利散射相函数，参考书籍：《大气辐射导论》 廖国南,P95*/
{
	return 0.75*(1+cos(x)*cos(x));
}

double HGspfun(double x,double g)
/*HG散射相函数，参考文献：大气气溶胶粒子散射相函数的数值计算[J],红外与激光工程*/
{
	return 1/(4*PI)*(1-g*g)/pow(1+g*g-2*g*cos(x),1.5);
}

double IHGspfun(double x,double g)
/*改进的HG散射相函数，参考文献：大气气溶胶粒子散射相函数的数值计算[J],红外与激光工程*/
{
	return 1/(4*PI)*1.5*(1-g*g)*(1+cos(x)*cos(x))/(pow(1+g*g-2*g*cos(x),1.5)*(2+g*g));
}
#undef  PI