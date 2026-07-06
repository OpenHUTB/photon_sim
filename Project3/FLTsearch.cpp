#include "stdafx.h"
#include <cmath>

int FLTsearch(double data[],int N,double X)
// 在数组data[]中查找与X最接近的数，返回其对应的下标号
{
	int index;
	double small;
	small = fabs(X-data[0]);
	index = 0;
	for(int i=0;i<N;i++){
		if(fabs(X-data[i])<small){
			small = abs(X-data[i]);
			index = i;
		}
     }
	return index;
}
