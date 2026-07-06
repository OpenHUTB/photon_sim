#include "stdafx.h"
#include <cmath>

#define PI         3.1415926
#define  MPX       1024
#define  MPY       1024
#define  MaxDPZA   100								/* 下行路径上观测天顶角数目的最大值 */
#define  MaxDPAA   100								/* 下行路径上观测方位角数目的最大值 */

void ImgSim(double WL_c,double Sensor_Height,double SunZA,double ViewZA,double ViewAA,char *PSF_fname,int ImageXN,int ImageYN)
{
	/* wavelength range for Radiance Parameter for Source/Sensor/Downwelled Path Files 
	对于源、传感器、下行路径文件的辐射参数范围 */
	double   WL_Begin =  500.0;                    /* the begin of the wavelength [nm] */
	double   WL_End =  600.0;                      /* the end of the wavelength [nm] */
	double   WL_Inc = 10.0;                        /* the increment of the wavelength [nm] */

	double   AlbedoT_Begin = 0.05;                 /* begin of the target albedo 目标反照率的开始值 */
	double   AlbedoT_Inc = 0.6;                    /* increment of the targe albedo 目标反照率的增量 */
	double   AlbedoBG_Begin = 0.05;                /* begin of the background albedo 背景反照率的开始值 */
	double   AlbedoBG_Inc = 0.3;                   /* increment of the background albedo 背景反照率的增量 */
	int      NT = 2;                               /* number of the target albedo 目标反照率的数目 */
	int      NBG = 3;                              /* number of the background albedo 背景反照率的数目 */

	double   SPSH_Begin = 100.0;                   /* begin of the sensor height in Sensor Path [km] 在传感器路径上传感器高度的开始值 */
	double   SPSH_End = 100.0;                     /* end of the sensor height in Sensor Path [km] 结束值 */
	double   SPSH_Inc = -1.0;                      /* increment of the sensor height in Sensor Path [km] 增量 */

	double   SPZA_Begin = 100.0;                   /* begin of the view zenith angle in Sensor Path [deg] 在传感器路径上观测天顶角的起始值 */
	double   SPZA_End = 180.0;                     /* end of the view zenith angle in Sensor Path [deg] 结束值 */
	double   SPZA_Inc = 5.0;                       /* increment of the view zenith angle in Sensor Path [deg] 增量*/
	double   SPAA_Begin = -180.0;                  /* begin of the view azimuth angle in Sensor Path [deg] 传感器路径上观测方位角的起始值 */
	double   SPAA_End = 170.0;                     /* end of the view azimuth angle in Sensor Path [deg] 结束值 */
	double   SPAA_Inc = 10.0;                      /* increment of the view azimuth angle in Sensor Path [deg] 增量 */
	int      SPZA_Nc, SPAA_Nc, SPSH_Nc;
	int      N_SPZA, N_SPAA, N_SPSH;               /* number of the view zenith angle, azimuth angle and sensor height in Sensor Path 在传感器路径上观测天顶角、方位角、传感器高度的数目 */

	double   DPZA_Begin = 0.0;                     /* begin of the view zenith angle in Downwelled Path 在下行路径上观测天顶角的起始值 */
	double   DPZA_End = 84.0;                      /* end of the view zenith angle in Downwelled Path 结束值 */
	double   DPZA_Inc = 12.0;                      /* increment of the view zenith angle in Downwelled Path 增量 */
	double   DPAA_Begin = -180.0;                  /* begin of the view azimuth angle in Downwelled Path 在下行路径上观测方位角的起始值 */
	double   DPAA_End = 160.0;                     /* end of the view azimuth angle in Downwelled Path 结束值 */
	double   DPAA_Inc = 20.0;                      /* increment of the view azimuth angle in Downwelled Path 增量 */
	int      N_DPZA, N_DPAA;                       /* number of the view zenith angle, azimuth angle in Downwelled Path 下行路径上观测天顶角、方位角的数目 */
	
	double   AA;
	double   *ZA = new double[MaxDPZA];				/* */
	double   ** SR = new double *[MaxDPZA];        /* array of downwelled radiance from different view zenith angle and azimuth angle 从不同的观测天顶角和方位角的下行辐射的数组 */

	ifstream PSFFile,SourcePF,SensorPF,DownPF;		/* 点扩散输入文件、源输入文件、传感器输入文件、下行路径输入文件 */
    ofstream DownRad,RadF;							/* 下行辐射输出文件 RadF ???*/	
	int      N_WL;
	double   Albedo_T, Albedo_BG;					/* 目标反照率、背景反照率 */
	double   TotalRad,SolarRad,SkyRad,PathRad,AdjRad,BeamEn,TarRad;	/* 整个辐射、太阳辐射A、天空辐射B、路径辐射C+F、临近辐射I、入射光？、目标辐射 */
	int      ** PSF = new int *[MPX];				/*点扩散函数*/
	double   ** DPSF = new double *[MPX];			/**/
	double   ** RadPL = new double *[MPX];			/* ???*/
	double   ** Alb = new double *[MPX];			/* 反照率 */
	double   X,Y;                     /* [m] */
	int      i,j,m,n;
	int      TXI, TYI;                             /* index of target pixel in the image array 在图像数组中目标像素的的索引 */
	double   WL,TransUp,TransDown,SolarE;			/* 波长、上行传输、下行传输、太阳能量 */
	long     SourcePF_offset, SensorPF_offset, DownPF_offset1, DownPF_offset2;			/* 源参数文件偏移、传感器参数文件偏移、下行参数文件偏移1、2*/
	char     Rad_fname[50]={0};						/* 辐射文件名 */

	for(i=0;i<MaxDPZA;i++) {
		SR[i] = new double[MaxDPAA];
	}
	for(i=0;i<MPX;i++) {
		PSF[i] = new int[MPY];
		DPSF[i] = new double[MPY];
		RadPL[i] = new double[MPY];
		Alb[i] = new double[MPY];
	}

	PSFFile.open(PSF_fname);						/* 点扩散函数文件 */
	if (!PSFFile) {
		cerr << "error: unable to open " << PSF_fname << " file" << endl;
		exit(1);
	}
	
	SourcePF.open("Radiance Parameter for Source Path");				/* 源参数文件 */
	if (!SourcePF) {
		cerr << "error: unable to open input file" << endl;
		exit(1);
	}

	SensorPF.open("Radiance Parameter for Sensor Path");				/* 传感器参数文件 */
	if (!SensorPF) {
		cerr << "error: unable to open input file" << endl;
		exit(1);
	}

	DownPF.open("Radiance Parameter for Downwelled Path");				/* 下行参数文件 */
	if (!DownPF) {
		cerr << "error: unable to open input file" << endl;
		exit(1);
	}

    DownRad.open("Downwelled Radiance.dat",ios_base::app);				/* 打开下行辐射文件 */
	    if (!DownRad) {
		    cerr << " error: unable to open input file" << endl;
		    exit(1);
	    }

	for(i=0;i<ImageXN;i++) {
		for(j=0;j<ImageYN;j++) {
			PSFFile >> m >> n >> PSF[i][j] >> DPSF[i][j];
		}
	}

	TXI = ImageXN/2;
	TYI = ImageYN/2;

/* Open the radiance file and read the radiance parameters for specific wavelength 对特定的波长，打开辐射文件和辐射参数 */		
	/* read radiance parameters for the source path 对源路径，读取辐射参数 */
	SourcePF_offset = (WL_c-WL_Begin)/WL_Inc*(3*8+7);
	SourcePF.seekg(SourcePF_offset,ios::beg);
	SourcePF >> WL >> SolarE >> TransDown;
	// 注意 太阳方向和表面法线方向之间的夹角
	SolarRad = SolarE*TransDown*cos(SunZA/180.0*PI);  //pay attention to the angle between the direction of the sun and the direction of the surface normal  
		
	/* read radiance parameters for the sensor path 读取传感器路径上的辐射参数 */
	//Attention:according to given zenith angle and azimuth angle,find path radiance and upwelled transmission in the closest direction, rather than interpolating method 
	N_SPZA = int(fabs(SPZA_End - SPZA_Begin)/fabs(SPZA_Inc) + 1 );
	N_SPAA = int(fabs(SPAA_End - SPAA_Begin)/fabs(SPAA_Inc) + 1 );
	N_SPSH = int(fabs(SPSH_End - SPSH_Begin)/fabs(SPSH_Inc) + 1 );
	N_WL = int(fabs(WL_End - WL_Begin)/fabs(WL_Inc) + 1);
	SPSH_Nc = int((Sensor_Height - SPSH_Begin)/SPSH_Inc);
	if(fabs(Sensor_Height-(SPSH_Nc*SPSH_Inc+SPSH_Begin)) > fabs(Sensor_Height-((SPSH_Nc+1)*SPSH_Inc+SPSH_Begin))) SPSH_Nc += 1;
	SPZA_Nc = int((ViewZA - SPZA_Begin)/SPZA_Inc);
	if(fabs(ViewZA-(SPZA_Nc*SPZA_Inc+SPZA_Begin)) > fabs(ViewZA-((SPZA_Nc+1)*SPZA_Inc+SPZA_Begin))) SPZA_Nc += 1;
	SPAA_Nc = int((ViewAA - SPAA_Begin)/SPAA_Inc);
	if(fabs(ViewAA-(SPAA_Nc*SPAA_Inc+SPAA_Begin)) > fabs(ViewAA-((SPAA_Nc+1)*SPAA_Inc+SPAA_Begin))) SPAA_Nc += 1;
	SensorPF_offset = (SPSH_Nc*N_SPZA*N_SPAA + SPZA_Nc*N_SPAA + SPAA_Nc)*(3*8+8+N_WL*(3*8+7))+(3*8+8)+(WL_c-WL_Begin)/WL_Inc*(3*8+7);

	SensorPF.seekg(SensorPF_offset,ios::beg);
	SensorPF >> WL >> PathRad >> TransUp;

	/* read radiance parameters for the downwelled path and calculate the sky radiance by integral */
	DownPF_offset1 = (WL_c-WL_Begin)/WL_Inc*(2*8+8);
	DownPF_offset2 = (WL_End-WL_c)/WL_Inc*(2*8+8);
    N_DPZA = int(fabs(DPZA_End - DPZA_Begin)/fabs(DPZA_Inc) + 1);
	N_DPAA = int(fabs(DPAA_End - DPAA_Begin)/fabs(DPAA_Inc) + 1);
    DownRad << WL_c << endl;
	DownPF.seekg(0,ios_base::beg);
	for (i=0;i<N_DPZA;i++) {
		for(j=0;j<N_DPAA;j++) {
			DownPF >> ZA[i] >> AA;
			DownPF.seekg(DownPF_offset1,ios::cur);
			DownPF >> WL >> SR[i][j];
			DownPF.seekg(DownPF_offset2,ios::cur);
			DownRad.width(5);
			DownRad << ZA[i];
			DownRad.width(5);
			DownRad << AA;
			DownRad.width(10);
			DownRad << SR[i][j] << endl;
		}
	}

    /* calculate the sky radiance i.e.downwelled radiance */
	SkyRad = 0.0;
	for (i=0;i<N_DPZA;i++) {
		for(j=0;j<N_DPAA;j++) {
			SkyRad += SR[i][j]*cos(ZA[i]/180.0*PI)*sin(ZA[i]/180.0*PI)*(DPZA_Inc/180.0*PI)*(DPAA_Inc/180.0*PI);
		}
	}

/* calculate the total radiance at the sensor according to different albedo for target and background */
	for (m=0;m<NT;m++) {
		for (n=0;n<NBG;n++) {
            /* create the albedo array */
			Albedo_T = AlbedoT_Begin + m*AlbedoT_Inc;
			Albedo_BG = AlbedoBG_Begin + n*AlbedoBG_Inc;
			cout << "   Target Albedo = " << Albedo_T;
			cout << ",  Background Albedo = " << Albedo_BG << endl;
			for (i=0;i<ImageXN;i++) {
				for (j=0;j<ImageYN;j++) {
					Alb[i][j] = Albedo_BG;
				}
			}
			if (ImageXN%2==0) {
				if (ImageYN%2==0) {
					Alb[TXI][TYI] = Albedo_T;
					Alb[TXI-1][TYI] = Albedo_T;
					Alb[TXI][TYI-1] = Albedo_T;
					Alb[TXI-1][TYI-1] = Albedo_T;
				}
				else {
					Alb[TXI-1][TYI] = Albedo_T;
					Alb[TXI][TYI] = Albedo_T;
				}
			}
			else {
				if (ImageYN%2==0) {
					Alb[TXI][TYI-1] = Albedo_T;
					Alb[TXI][TYI] = Albedo_T;

				}
				else {
					Alb[TXI][TYI] = Albedo_T;
				}
			}

            /* create the surface leaving radiance array */
			for (i=0;i<ImageXN;i++) {
				for (j=0;j<ImageYN;j++) {
					RadPL[i][j] = (SolarRad + SkyRad) * Alb[i][j] / PI;
				}
			}
            /* calculate the beam enhancement and the radiance reflected by the target area */
			if (ImageXN%2==0) {
				if (ImageYN%2==0) {
					BeamEn = RadPL[TXI][TYI]*DPSF[TXI][TYI] + RadPL[TXI-1][TYI]*DPSF[TXI-1][TYI] + RadPL[TXI][TYI-1]*DPSF[TXI][TYI-1] + RadPL[TXI-1][TYI-1]*DPSF[TXI-1][TYI-1];
					TarRad = TransUp*(RadPL[TXI][TYI]+RadPL[TXI-1][TYI]+RadPL[TXI][TYI-1]+RadPL[TXI-1][TYI-1]);
				}
				else {
					BeamEn = RadPL[TXI-1][TYI]*DPSF[TXI-1][TYI] + RadPL[TXI][TYI]*DPSF[TXI][TYI];
					TarRad = TransUp*(RadPL[TXI-1][TYI]+RadPL[TXI][TYI]);
				}
			}
			else {
				if (ImageYN%2==0) {
					BeamEn = RadPL[TXI][TYI-1]*DPSF[TXI][TYI-1] + RadPL[TXI][TYI]*DPSF[TXI][TYI];
					TarRad = TransUp*(RadPL[TXI][TYI-1]+RadPL[TXI][TYI]);
				}
				else {
					BeamEn = RadPL[TXI][TYI]*DPSF[TXI][TYI];
					TarRad = TransUp*RadPL[TXI][TYI];
				}
			}
            /* calculate the cross radiance */
			AdjRad = 0.0;
			BeamEn = 0.0;
			for (i=0;i<ImageXN;i++) {
				for (j=0;j<ImageYN;j++) {
					AdjRad += RadPL[i][j] * DPSF[i][j];
				}
			}
			AdjRad -= BeamEn;

			/* calculate the total radiance at the sensor */
			TotalRad = TarRad + BeamEn + AdjRad + PathRad;
			cout << "     Total Radiance at the sensor is: " << TotalRad << endl;
			cout << "     Beam enhancement is: " << BeamEn << endl;
			cout << "     Radiance caused by Adjacent effect is: " << AdjRad << endl;
			cout << "     Path radiance is: " << PathRad << endl;

            /*open files to write calculation results */
			sprintf(Rad_fname, "Radiance_T%02dBK%02d.dat",int(Albedo_T*100),int(Albedo_BG*100));
			RadF.open(Rad_fname,ios_base::app);
			if (!RadF) {
				cerr << "error: unable to open input file" << endl;
				exit(1);
			}
			if(WL_c == WL_Begin) {
				RadF <<"Target Albedo = " << Albedo_T << "    Background Albedo = " << Albedo_BG << endl;
				RadF <<"WL" << "             TotalRad" << "       TargetRad" << "      BeamEn" << "         AdjRad" << "         PathRad" << "        A/T(%)"<< endl;
			}
			RadF.width(8);
			RadF << scientific << WL;
			RadF.width(15);
			RadF << right << scientific << TotalRad;
			RadF.width(15);
			RadF << scientific << TarRad;
			RadF.width(15);
			RadF << scientific << BeamEn;
			RadF.width(15);
			RadF << scientific << AdjRad;
			RadF.width(15);
			RadF << scientific << PathRad;
			RadF << "  ";
			RadF << left  << fixed << AdjRad/TotalRad*100 << endl;
			RadF.close();
		}
	}

    DownRad.close();
	DownPF.close();
	SensorPF.close();
	SourcePF.close();
	PSFFile.close();

	for (i=0;i<MaxDPZA;i++) {
		delete SR[i];
		SR[i] = NULL;
	}
	delete[] SR;
	SR = NULL;

	for(i=0;i<MPX;i++) {
		delete PSF[i];
		PSF[i] = NULL;
		delete DPSF[i];
		DPSF[i] = NULL;
		delete RadPL[i];
		RadPL[i] = NULL;
		delete Alb[i];
		Alb[i] = NULL;
	}
	delete[] PSF;
	PSF = NULL;
	delete[] DPSF;
	DPSF = NULL;
	delete[] RadPL;
	RadPL = NULL;
	delete[] Alb;
	Alb = NULL;
}
#undef   MPX
#undef   MPY
#undef   MaxDPZA
#undef   MaxDPAA
#undef   PI          
