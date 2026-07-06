#include "stdafx.h"
#include <cmath>

#define	 PI       3.1415926535
#define  MaxNP    10000000
#define  MPX      1024
#define  MPY      1024

/* 点扩散函数统计：图像X轴方向的像素数目；图像Y轴方向的像素数目; 网格的宽度 1米 [m]*/
void PSFStatistic(int ImageXN,int ImageYN,double GWidth,int Nphoton,char *Pos_fname,char *PSF_fname)
{
	ifstream PosFile;
	ofstream PSFFile;

	double   *PosX, *PosY;                           /* [km] */
	double   X, Y;
	double   GridWidth, GridLength;
	int      SumM, SumA;
	int      i, j;
	int      SumP;
	double   Width, Length;

	int ** PSF = new int *[MPX];
	double ** DPSF = new double *[MPX];
	for (i=0;i<MPX;i++) {
		PSF[i] = new int[MPY];
		DPSF[i] = new double[MPY];
	}
	
	PosX = new double[MaxNP];
	PosY = new double[MaxNP];
	GridWidth = GWidth;
	GridLength = GWidth;

	PosFile.open(Pos_fname);
	if (!PosFile) {
		cerr << "error: unable to open " << Pos_fname << " file" << endl;
		exit(1);
	}
	for(i=0;i<Nphoton;i++) {
		PosFile >> PosX[i] >> PosY[i] >> SumM >> SumA;
	}
    
	Width=ImageXN*GridWidth;		// 图像X轴方向的像素数目*网格宽度
	Length=ImageYN*GridLength;

	for(i=0;i<ImageXN;i++) {
		for(j=0;j<ImageYN;j++) {
			PSF[i][j] = 0;
			DPSF[i][j] = 0.0;
		}
	}

    for(i=0;i<Nphoton;i++) {
		X=1000.0*PosX[i];              //turn km to m
		Y=1000.0*PosY[i];
		X += Width/2;
		Y += Length/2;
		if(0<=X && 0<=Y) {
			if(X<=Width && Y<=Length) {
				if(X==Width) {
					if(Y==Length) {
						PSF[ImageXN-1][ImageYN-1] += 1;
					}
					else {
						PSF[ImageXN-1][int(Y/GridLength)] += 1;
					}
				}
				else {
					if(Y==Length) {
						PSF[int(X/GridWidth)][ImageYN-1] += 1;
					}
					else {
							PSF[int(X/GridWidth)][int(Y/GridLength)] += 1;
					}
				}
			}
		}
	}

	SumP = 0;
	for(i=0;i<ImageXN;i++) {
		for(j=0;j<ImageYN;j++) {
            DPSF[i][j] = 1.0*PSF[i][j]/Nphoton;		// 光子数目占所有的百分比
			SumP += PSF[i][j];		// 光子总数目统计
		}
	}
	cout << "   Summary of the photons in the image: " << SumP << endl;

	PSFFile.open(PSF_fname);
	if (!PSFFile) {
		cerr << "error: unable to open " << PSF_fname << " file" << endl;
		exit(1);
	}

	for(i=0;i<ImageXN;i++) {
		for(j=0;j<ImageYN;j++) {
			PSFFile.width(5);
			PSFFile << i;
			PSFFile.width(5);
			PSFFile << j;
			PSFFile.width(10);
			PSFFile << PSF[i][j];	// (i,j)位置上有多少个点（点扩散而成）
			PSFFile.width(20);
    		PSFFile << DPSF[i][j] << endl;			// (i,j)位置上的光子数占全部数目的百分比		
		}
	}

	PSFFile.close();
	PosFile.close();

	delete[] PosX;
	delete[] PosY;
	
	for(i=0;i<MPX;i++) {
		delete PSF[i];
		PSF[i] = NULL;
		delete DPSF[i];
		DPSF[i] = NULL;
	}
	delete[] PSF;
	PSF = NULL;
	delete[] DPSF;
    DPSF = NULL;
}
#undef   MPX
#undef   MPY
#undef	 PI        
#undef   MaxNP
