#include "stdafx.h"
#include <cmath>
#include <time.h>
#define  MaxH        100
#define  MaxWL       1000
#define  MaxNP       10000000
#define  MaxSPF      180		/* 最大散射相函数 */
#define	 PI          3.1415926535
#define  ALIVE       1   		/* if photon not yet terminated 光子仍或者 alive*/
#define  DEAD        0    		/* if photon is to be terminated 光子死了 dead */
#define  THRESHOLD   0.01		/* 门限值 */
#define  ONE_MINUS_COSZERO 1.0E-12
     /* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
     /* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */
#define SIGN(x)           ((x)>=0 ? 1:-1)

/**
WL: 目前的波长；N_H2： 垂直高度；遥感器高度；观测天顶角；观测方位角；光子数目；Pos文件名；网格宽度；sa: 散射角度；Mnspf: 归一化的分子散射相函数； Anspf: 归一化的气溶胶散射相函数
*/
void MonteCarlo(double WL,int N_H2,double Sensor_Height,double PT_ZA,double PT_AA,int Nphoton,char *Pos_fname,double GridWidth,double sa[],double Mnspf[],double Anspf[], char* OpticalDepthFileName)
{
	/* wavelength range for Optical Depth File */
	double   WL_Begin =  500.0;              /* the begin of the wavelength [nm] */
	double   WL_End =  600.0;                /* the end of the wavelength [nm] */
	double   WL_Inc = 10.0;                  /* the increment of the wavelength [nm] */

	ifstream  ODSFile;                       /* Original Data file(altitude,wavelength,optical depth) 原始数据（高度、气溶胶光学厚度、分子光学厚度、总光学厚度）*/
	ofstream  PosFile;
	ofstream  ODFile;

	char      OD_fname[50] = {0};

	double    *Alt;                          /* array of the vertical altitude [km] 垂直高度的数组 */
	double    *Mol;                          /* array of the molucule optical depth along the vertical altitude 沿着垂直高度方向 分子光学厚度 的数组*/
	double    *Aer;                          /* array of the aerosol optical depth along the vertical altitude 沿着垂直高度方向 气溶胶光学厚度 的数组*/
	double    *Tol;                          /* array of the total optical depth along the vertical altitude 沿着垂直高度方向 整个光学厚度 的数组*/

	double    *PosX, *PosY;                  /* photon position in x and y directions after propagation [km] 传播之后，光子在在x轴和y轴方向上的位置 数组*/
	int       *SumM;                         /* collision times of every photon with molucule(molecule) 每一个光子和分子碰撞的次数 */
	int       *SumA;                         /* collision times of every photon with aerosol 每一个光子和气溶胶碰撞的次数*/
	int       i_photon;                      /* current photon 目前的光子 */
	int       CTM;                           /* collision times of photon with molucule 光子和分子碰撞的次数*/
	int       CTA;                           /* collision times of ephoton with aerosol 光子和气溶胶碰撞的次数 */
	double    W;                             /* photon weight 光子的权重 */
	short     photon_status;                 /* flag = ALIVE=1 or DEAD=0 光子的状态，或生或死*/
	double    x, y, z;                       /* photon position [m] 光子的位置 */
	double    ux, uy, uz;                    /* photon trajectory as cosines 光子的飞行轨迹 ？？？？？？*/
	double    uxx, uyy, uzz;                 /* temporary photon trajectory values 临时光子轨道的值 ？？？？？？？？？？*/
	double	  ds;                            /* step size in photon trajectory direction 在光子的轨迹方向上移动的步子的大小 */
	double    dz;                            /* step size in vertical direction 在垂直方向上移动步子的大小 */
	double    ol;                            /* optical length 光程？？？？？？？？？？？*/
	double    r_theta;                       /* random number decide photon launch direction 随机数决定光子发射的方向 */
	double    r_psi;                         /* random number decide photon launch direction */
	double    r_ol;                          /* random number decide optical length 随机数决定光程 */
	double    r_cp;                          /* random number decide the collision proterty 决定碰撞的??????*/
	double    r_sa;                          /* random number decide the scatter angle 随机数决定散射的角度 */
	double    r_raa;                         /* random number decide the scatter relative azimuth angle 随机数决定散射的相对方位角 */
	double    cp;                            /* collision probality 碰撞的概率 */
	double    psi;                           /* azimuth angle 方位角 */
    double	  costheta;                      /* cos(theta) */
	double    sintheta;                      /* sin(theta) */
	double    cospsi;                        /* cos(psi) 方位角的余弦 */
	double    sinpsi;                        /* sin(psi) */
	int       index, index1, index2;
	double    trans;                         /* transmittance between two positions in photon trajectory (slant path, not vertical path)在光子轨道上（倾斜的）两点之间的透光率 */
	double    aod;                           /* aerosol optical depth 气溶胶光学厚度 */
	double    mod;                           /* molecule optical depth 分子光学厚度 */
	double    tod;                           /* total optical depth 整个光学厚度 */
	double    c_sa;                          /* current scatter angle after collision [rad] 碰撞之后当前的散射角度 */
	double    raa;                           /* scatter relative azimuth angle after collision [rad] 碰撞之后散射的相对方位角 */
	double    GridWidth2, half_angle;
	double    MaxZA, MinZA, MaxAA, MinAA;
	double    dtod, tod1, tod2;
	double    tempd;
	long      offset1, offset2;
	int       i;
	double    Angle[4];
	double    NA = 4;
    
	// DECLARE FUNCTION 声明外部函数
	extern int FLTsearch(double [],int ,double);
	extern double RandomGen(char , long , long *);     // Random number generator
	extern void EmitAngle(double [], int , double , double , double , double);

	// Obtain system current time as the seed of the random generator 获得系统当前时间作为随机数产生器的种子
	time_t rawtime;
	tm *pTm;
	time(&rawtime);
	long RSeed;
	pTm = localtime(&rawtime);
	RSeed = (pTm->tm_year+1900)+(pTm->tm_mon+1)+(pTm->tm_mday)+(pTm->tm_hour)+(pTm->tm_min)+(pTm->tm_sec);	/* 随机种子 */

	PT_ZA = 180.0 - PT_ZA;// MODTRAN中的 天顶角 与 大倾角观测时的观测角 之间的转换
	PT_ZA *= (PI/180.0);                                             // turn [degree] to [radian] 将[度]转换为[弧]
	PT_AA *= (PI/180.0);                                             // turn [degree] to [radian]
	Sensor_Height *= 1000.0;                                         // turn [km] to [m]
	GridWidth2 = GridWidth / 2.0;								     // 观测长度范围的一半 / 高度 = 半角(half_angle)的正切值
	half_angle = atan( GridWidth2/Sensor_Height);                    // [ radian ]

	Alt = new double[MaxH];		// 垂直高度
	Mol = new double[MaxWL];	// 沿着垂直高度方向 分子光学厚度
	Aer = new double[MaxWL];	// 沿着垂直高度方向 气溶胶光学厚度
	Tol = new double[MaxWL];	// 沿着垂直高度方向 整个光学厚度

	PosX = new double[MaxNP];	// 传播之后，光子在在x轴和y轴方向上的位置
	PosY = new double[MaxNP];
	SumM = new int[MaxNP];		// 每一个光子和分子碰撞的次数
	SumA = new int[MaxNP];		// 每一个光子和气溶胶碰撞的次数


	for (i=0;i<MaxH;i++) {		// 垂直高度初始化
		Alt[i] = 0.0;
	}
	for (i=0;i<MaxWL;i++) {		// 光学厚度初始化
		Mol[i] = 0.0;
		Aer[i] = 0.0;
		Tol[i] = 0.0;
	}

    /*****  read altitude, molecule optical depth, aerosol optical depth, total optical depth from "Optical Depth" file *****/
	ODSFile.open(OpticalDepthFileName);
	//ODSFile.open("Optical Depth"); // 从"Optical Depth"中读取 高度(500-600)、分子光学厚度、气溶胶光学厚度、全部光学厚度(分子光学厚度+气溶胶光学厚度) ODS 读文件
	if (!ODSFile) {
		cerr << "error: unable to open input file" << endl;
		exit(1);
	}
	
	sprintf(OD_fname, "OD%04d.dat", int(WL));	// 光学厚度文件名：OD0500.dat ~ OD0600.dat
	ODFile.open(OD_fname);						// OD 写文件
	if (!ODFile) {
		cerr << "error: unable to open " << OD_fname << " file" << endl;
		exit(1);
	}

	i = 0;
	do {
		ODSFile >> Alt[i];					 // 读入 垂直高度
		Alt[i] *= 1000.0;                    // turn [km] to [m]
		ODFile.width(11);
		ODFile << Alt[i];					 // 将[km]转换成[m]后写入 垂直高度

		offset1 = (WL - WL_Begin) / WL_Inc * (4*8+14);	// 指针向后移动14个位置; (WL - WL_Begin) / WL_Inc: 第几行 510(第一行), 520(第二行)
		ODSFile.seekg(offset1, ios::cur);
		ODSFile >> WL >> Mol[i] >> Aer[i] >> Tol[i];
		ODFile.width(11);
		ODFile << Mol[i];
		ODFile.width(11);
		ODFile << Aer[i];
		ODFile.width(11);
		ODFile << Tol[i] << endl;            //要求给出 垂直光学厚度，建立垂直光学厚度与垂直高度的一一对应关系

		offset2 = (WL_End - WL) / WL_Inc * (4*8+14);
		ODSFile.seekg(offset2, ios::cur);

		i += 1;
	} while (i < N_H2);						// 将垂直高度分成的段数 (29.5-0.0)/(10+1)

	/**** INITIALIZATIONS ****/
	RandomGen(0, RSeed, NULL);               // Initializes the seed for the random number generator 用随机数产生器初始化随机种子

	/**** RUN
	Do Monte Carlo simulation
	****/
	i_photon = 0;							 // 第i个光子
	do {
		i_photon += 1;	                     // increment photon count
		CTM = 0;							 // 光子与分子的碰撞次数
		CTA = 0;							 // 光子与气溶胶的碰撞次数
//		if (10*i_photon==Nphoton) {
//			cout << "10.00% finished!" << endl;
//		}
//		else if (2*i_photon==Nphoton) {
//			cout << "50.00% finished!" << endl;
//		}
//		else if (i_photon==Nphoton) {
//			cout << "100.00% finished!" << endl;
//		}

		W = 1.0;                             // set photon weight to one 设置光子的权重为1
		photon_status = ALIVE;               // Launch an ALIVE photon 发射的是一个活着的光子

		// Set position of photon emission.	// 设置光子发射的位置
		x = Sensor_Height*tan(PT_ZA)*cos(PT_AA);        // [ m ]         
		y = Sensor_Height*tan(PT_ZA)*sin(PT_AA);        // [ m ]
		z = Sensor_Height;                              // [ m ]

		if ( fabs(PT_AA-0.5*PI)<ONE_MINUS_COSZERO || fabs(PT_AA-1.5*PI)<ONE_MINUS_COSZERO ) {
			x = 0.0;
		}
		if (PT_AA == PI || PT_AA == 2*PI) {
			y = 0.0;
		}

		// generate photon trajectory from source, the Instantaneous Field of View(IFOV) was considered 从源产生光子轨道，并考虑上瞬时视场
		// 计算发射光子天顶角和相对方位角的最大值和最小值
		EmitAngle(Angle,NA, x, y, z, GridWidth);      // calculate the maximum and the minimum of zenith angle and relative azimuth angle of the emitting photon  

		// determine the zenith angle of the emitting photon 决定发射光子的天顶角
		MaxZA = Angle[0];
		MinZA = Angle[1];
		r_theta = MinZA + RandomGen(1, 0, NULL)*( MaxZA - MinZA );
		costheta = cos(r_theta);  
		sintheta = sqrt(1.0 - costheta*costheta);	// sintheta is always positive

		// determine the relative azimuth angle of the emitting photon 决定发射光子的相对方位角
		MaxAA = Angle[2];
		MinAA = Angle[3];
		r_psi = RandomGen(1, 0, NULL);
		psi = MinAA + r_psi*(MaxAA-MinAA);

		ux = sintheta*cos(psi);
		uy = sintheta*sin(psi);
		uz = costheta;

        //calculate the photon position without atmosphere effect 没有大气的作用时计算光子的位置
/*		ds = Sensor_Height / uz;
		x += ds * ux;
		y += ds * uy;
		z -= ds * uz;
		PosX[i_photon-1] = x;
		PosY[i_photon-1] = y;
*/
		
		/**** MOVE_COLLIDE_ 移动相碰撞 ****/
		do {
			// MOVE
			r_ol = RandomGen(1, 0, NULL);
			ol = -log(1-r_ol);                         // generate random optical distance(i.e. optical depth in photon trajectory direction) 产生随机光学距离（比如：在光子轨迹方向上的光学厚度）
			dtod = ol * costheta;                      // calculate the increment of the optical depth in perpendicular direction 计算在垂直方向上的光学厚度的增量
			index = FLTsearch(Alt,N_H2,z);			   // 在数组data[](Alt[])中查找与X(z)最接近的数，返回其对应的下标号index
			tod1 = Tol[index];                         // look up last position's optical depth 查询最后位置的光学厚度 
			tod2 = tod1 + dtod;                        // calculate current position's optical depth 计算当前位置的光学厚度
            index1 = FLTsearch(Tol,N_H2,tod2);         // look up the height of the atmosphere layer of current position 查找当前位置大气层的高度
			dz = z - Alt[index1];                      // calculate altitude increment in vertical direction 计算在垂直方向上的高度增量
			if (1 - fabs(uz) <= ONE_MINUS_COSZERO)     // close to perpendicular. 接近于垂直
				ds = dz;							   // ds:在光子的轨迹方向上移动的步子的大小; dz:在垂直方向上移动步子的大小
			else 
				ds = dz / costheta;                    
			x += ds * ux;                              // Update positions. 更新光子的位置
			y += ds * uy;
			z -= ds * uz;
			if (z <= 0.0) {                           // photon arrives the ground, photon dead 光子到达地面就死了
                photon_status = DEAD;
				PosX[i_photon-1] = x;				  // PosX, PosY: 传播之后，光子在在x轴和y轴方向上的位置
				PosY[i_photon-1] = y;
			} else {
				trans = exp(-1.0*ol);                  // trans: 在光子轨道上（倾斜的）两点之间的透光率
				W *= trans;                            // photon weight after moving 移动之后光子的权重(乘以透光率)
				if (W < THRESHOLD) {                   // if photon weight is less than threshold, photon dead 如果光子的权重小于门限值，则光子死
					photon_status = DEAD;
					PosX[i_photon-1] = x;
					PosY[i_photon-1] = y;
				} else {                               // if photon weight is larger than threshold, photon alive, photon will collide with atmosphere particles 
					// COLLIDE 如果光子的权重值大于门限值，则光子是活的，光子就会和大气粒子发生碰撞
					r_cp = RandomGen(1,0,NULL);
					mod = Mol[index1];					// 分子光学厚度
					aod = Aer[index1];					// 气溶胶光学厚度
					cp = mod / ( mod + aod );					      // collision property 碰撞的特性
					r_sa = RandomGen(1,0,NULL);						  // 随机数决定散射的角度
					if (0 <= r_cp && r_cp <= cp) {                    // photon will collide with molecule 	如果随机生成的r_cp值小于碰撞特性(cp)，光子将会和分子发生碰撞				
						//if(typespf!='R'){
						//	typespf = 'R';
						//}		
						index2 = FLTsearch(Mnspf,MaxSPF,r_sa);
						CTM += 1;
					} else if (cp < r_cp && r_cp <= 1) {	          // photon will collide with aerosol 	若果随机生成的r_cp值大于碰撞特性(cp)，光子将会和气溶胶发生碰撞		
						//if(Atypespf=='H'||Atypespf=='I'){
			            //  typespf = Atypespf;
						//}
						index2 = FLTsearch(Anspf,MaxSPF,r_sa);		  // 搜索Anspf中最接近随机产生的散射角度(r_sa)的下标
						CTA += 1;
					}
					c_sa = sa[index2];                                // obtain scatter angle after collision 获得碰撞后的散射角度
					r_raa = RandomGen(1,0,NULL);
					raa = 2*PI*r_raa;                                 // obtain scatter relative azimuth angle after collision 获得碰撞后散射的相对方位角

					costheta = cos(c_sa);
					sintheta = sqrt(1.0 - costheta*costheta);         // sqrt() is faster than sin(). 
					cospsi = cos(raa);
					if (raa < PI)
						sinpsi = sqrt(1.0 - cospsi*cospsi);           // sqrt() is faster than sin(). 
					else
						sinpsi = -sqrt(1.0 - cospsi*cospsi);
					
					// New trajectory 参考文献：唐军武. 海洋光学特性模拟与遥感模型: [博士论文]. 北京：中国科学院遥感应用研究所, 1999.
					if (1 - fabs(uz) <= ONE_MINUS_COSZERO) {       // close to perpendicular. 接近于垂直
						uxx = sintheta * cospsi;
						uyy = sintheta * sinpsi;
						uzz = costheta * SIGN(uz);                 // SIGN() is faster than division.
					} else {					                   // usually use this option
						tempd = sqrt(1.0 - uz * uz);
						uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / tempd + ux * costheta;
						uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / tempd + uy * costheta;
						uzz = -sintheta * cospsi * tempd + uz * costheta;
					}
					// Update trajectory
					ux = uxx;
					uy = uyy;
					uz = uzz;
				}
			} 
		} while (photon_status == ALIVE);	// 直到光子死亡

		SumM[i_photon-1] = CTM;		// 光子与分子碰撞的次数
		SumA[i_photon-1] = CTA;		// 光子与气溶胶的碰撞次数
	} while (i_photon < Nphoton);	// 直到把所有光子都模拟完

	//Output the position of the photon 输出光子的位置
	PosFile.open(Pos_fname);
	if (!PosFile) {
		cerr << "error: unable to open " << Pos_fname << " file" << endl;
		exit(1);
	}
	PosFile.setf(ios_base::scientific, ios_base::floatfield);
	for(i=0;i<Nphoton;i++) {
		PosFile.width(15);		// 占15位
		PosFile << PosX[i]/1000.0;                     // Positions saved in unit of km 位置保存的单位是km; 传播之后(光子死亡)，光子在在x轴和y轴方向上的位置（统计光子的位置）
		PosFile.width(15);
		PosFile << PosY[i]/1000.0;                     // Positions saved in unit of km
		PosFile.width(8);
		PosFile << SumM[i];		// 每一个光子和分子碰撞的次数
		PosFile.width(8);
		PosFile << SumA[i] << endl;		// 每一个光子和气溶胶碰撞的次数
	}

	ODSFile.close();
	PosFile.close();	
	ODFile.close();

	delete[] Alt;
	delete[] Mol;
	delete[] Aer;
	delete[] Tol;
	delete[] PosX;
	delete[] PosY;
	delete[] SumM;
	delete[] SumA;

}
#undef   MaxH  
#undef   MaxWL 
#undef   MaxNP 
#undef   MaxSPF
#undef 	 PI    
#undef   ALIVE    
#undef   DEAD     
#undef   THRESHOLD 
#undef   ONE_MINUS_COSZERO