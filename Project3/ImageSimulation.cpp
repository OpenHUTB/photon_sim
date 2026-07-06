// ImageSimulation.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <cmath>
#define  MaxSPF  180								/* 最大的 散射相函数 */

int _tmain(int argc, _TCHAR* argv[])
{
	double   WL_Begin =  500.0;                    /* the begin of the wavelength 波长开始 [nm] */
	double   WL_End =  600.0;                      /* the end of the wavelength 波长结束 [nm] */
	double   WL_Inc = 10.0;                        /* the increment of the wavelength [nm] */

	double   Sensor_Height = 100.0;                /* the height of the sensor 传感器的高度 [km] */
	double   H2_Begin = 29.5;                      /* the begin of the vertical altitude 垂直高度的开始 ？？？？？？？？？？[km] */
	double   H2_End = 0.0;                         /* the end of the vertical altitude [km] */
	double   H2_Inc = 0.5;                         /* the increment of the vertical altitude [km] */

	double   SunZA = 0.0;                          /* sun zenith angle 太阳的天顶角 [deg] */

	double   ViewZA = 180.0;                       /* view zenith angle 观察的天顶角 [deg] */
	double   ViewAA = 0.0;                         /* view azimuth angle 观察的方位角 [deg] */ 

	/*double SunZAInc = 12.0;*/							/**/
	//double ViewZAInc = 12.0;						/**/
	double SunZAInc = 5.0;
	double ViewZAInc = 5.0;							// 观测天顶角每次改变的量减少为5度

	int      ImageXN = 1024;                       /* pixel numbers in x direction of the image 图像X轴方向的像素数目 */
	int      ImageYN = 1024;                       /* pixel numbers in y direction of the image */
    double   GridWidth = 1.0;                      /* 网格的宽度 1米 [m] */

	int      Nphoton = 1000000;                    /* number of photons */

	char     Mtypespf;                             /* type of the scatter phase function:'R' raylight scatter phase function for molecule;
										                                                 'H' HG scatter phase function for aerosol;
																				         'I" Improved HG scatter phase function for aerosol 
																						 分子的散射相位函数的类型：瑞丽散射、HG散射、增强的HG散射 ？？？？？？？*/
	char     Atypespf;                             /* type of the aerosol phase function:'H' HG scatter phase function for aerosol;
																				         'I" Improved HG scatter phase function for aerosol 
																						 气溶胶的散射相位函数类型：HG散射、增强的HG散射*/
	double   g;                                    /* asymetry factor to use with the Henyey-Greenstein aerosol phase function 
												   GreenStein气溶胶相位函数的非对称因子 ？？？？？？？？？？？？？？？*/

	double   WL;                                   /* current wavelength 目前的波长 [nm] */
	double   *sa;                                  /* array of the scatter angle 散射角度矩阵 [rad] */
	double   *Mnspf;                               /* array of the normalized scatter phase function for molecule 
												   归一化的分子散射相位函数 */
	double   *Anspf;                               /* array of the normalized scatter phase function for aerosol 
												   归一化的气溶胶散射相位函数 */
	int      N_WL, N_H2;							/* N_WL：将波长分成的段数； N_H2：将垂直高度所分成的段数 */
	char     Pos_fname[50] = {0};
	char     PSF_fname[50] = {0};
	int      i;

	Mtypespf = 'R';
	// cout << "Input aerosol phase function type('H' for HG; 'I' for improved HG): ";
	// cin >> Atypespf;
	Atypespf = 'I';
	// cout << "Input asymetry factor to use by HG or IHG phase function: ";
	// cin >> g;
	g = 0.9;

	N_WL = int( fabs(WL_Begin - WL_End) / fabs(WL_Inc) +1 );
	N_H2 = int( fabs(H2_Begin - H2_End) / fabs(H2_Inc) +1 );

    extern void Nspfun(char, int, double [], double [],double );
	extern void MonteCarlo(double,int,double,double,double,int,char *,double,double [],double [],double [], char*);
	extern void PSFStatistic(int,int,double,int,char *,char *);
	extern void ImgSim(double,double,double,double,double,char *,int,int);

	sa = new double[MaxSPF];
	Mnspf = new double[MaxSPF];
	Anspf = new double[MaxSPF];
	for (i=0;i<MaxSPF;i++) {
		sa[i] = 0.0;
		Mnspf[i] = 0.0;
		Anspf[i] = 0.0;
	}
	// sa: 散射角度； Mnspf: 归一化的散射相函数？？？
    Nspfun(Mtypespf,MaxSPF,sa,Mnspf,0.0);    // obtain normalized molecule phase function 获得分子归一化的散射相函数
	Nspfun(Atypespf,MaxSPF,sa,Anspf,g);      // obtain normalized aerosol phase function 气溶胶 g:气溶胶相位函数的非对称因子

	char* OpticalDepthFileName[4*6] = {"Optical Depth12", "Optical Depth14", "Optical Depth16", "Optical Depth19", 
									"Optical Depth22", "Optical Depth24", "Optical Depth26", "Optical Depth29", 
									"Optical Depth32", "Optical Depth34", "Optical Depth36", "Optical Depth39", 
									"Optical Depth42", "Optical Depth44", "Optical Depth46", "Optical Depth49", 
									"Optical Depth52", "Optical Depth54", "Optical Depth56", "Optical Depth59", 
									"Optical Depth62", "Optical Depth64", "Optical Depth66", "Optical Depth69" };
	//WL = WL_Begin;	// 当前计算波长

	cout << " !!Execute Monte Carlo Simulation..." << endl;
	MonteCarlo(WL,N_H2,Sensor_Height,(180-30),90,Nphoton,Pos_fname,GridWidth,sa,Mnspf,Anspf);
	//void MonteCarlo(double WL,int N_H2,double Sensor_Height,double PT_ZA,double PT_AA,int Nphoton,char *Pos_fname,double GridWidth,double sa[],double Mnspf[],double Anspf[], char* OpticalDepthFileName)
	//MonteCarlo(WL,N_H2,Sensor_Height,viewZA,ViewAA,Nphoton,Pos_fname,GridWidth,sa,Mnspf,Anspf);	//模特卡罗模拟 太阳天顶角没用？？？？
	sprintf(PSF_fname, "PSF_%03d_%03d.dat", (int)(180-30), (int)(90));	// 输出到的文件名：PSF0500.dat ~ PSF0600.dat

	cout << " !!Calculate the Point Spread Function..." << endl;
	PSFStatistic(ImageXN,ImageYN,GridWidth,Nphoton,Pos_fname,PSF_fname);	// 点扩散函数的统计
	
	delete[] sa;
	delete[] Mnspf;
	delete[] Anspf;

	system("pause");
	return 0;
}

