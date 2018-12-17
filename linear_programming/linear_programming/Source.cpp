#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <vector>
#include <algorithm>
#include <set>
#include <iostream>
using namespace std;
#define N (1001)
#define INFTY (1e30)
#define EXP (1e-10)
//INFIY 正无穷
//EXP	无穷小量,主要为了对付double 在判断时可能出现误差
double co_matrix[N][N] = { 0 };
double b_matrix[N]= { 0 };
double c_matrix[N] = { 0 };
set<int> BI;//存储基本列的下标
double z;//存储最小值
double theta[N] = { 0 };
double lambda[N] = { 0 };
int INITLIZESIMPLEX(int m,int n);
vector<double > CALCULATEX(int m,int n);
void PIVOT(int i,int j,int m,int n);//以co_matrix的aij为主元,执行一次高斯约旦消元
int find_min_index(double * vec,int len)
{
	double tmp = vec[0];
	int index = 0;
	for (int i = 0; i < len; i++)
	{
		if (tmp>vec[i])
		{
			tmp = vec[i], index = i;//选第1个最小的
		}
	}
	return index;
}
int find_max_index(double * vec, int len)
{
	double tmp = vec[0];
	int index = 0;
	for (int i = 0; i < len; i++)
	{
		if (tmp<vec[i])
		{
			tmp = vec[i], index = i;//选第1个最大的
		}
	}
	return index;
}
pair < vector<double>, double > SIMPLEX(int m, int n, int init_flag = 0)
{
	vector<double> x;
	int feasible = INITLIZESIMPLEX(m, n);//可能有解的时候
	while (feasible)
	{
		int min_ci = find_min_index(c_matrix, n);
		if (c_matrix[min_ci] >= 0)
			//所有的ci都大于0,说明目标函数已经收敛了 找到了最优值
		{
			x = CALCULATEX(m,n);
			break;
		}
		//针对第min_ci列 计算最小的非负theta
		for (int i = 0; i < m; i++)
		{
			if (co_matrix[i][min_ci]>0)
				//这样计算的theta才是非负的
			{
				theta[i] = b_matrix[i] / co_matrix[i][min_ci];
			}
			else
				//当co_matrix[i][min_ci]为0,或者为负的时候,都不用去考虑
			{
				theta[i] = INFINITY;
			}
		}
		int max_thetai = find_min_index(theta, m);
		if (abs(theta[max_thetai] - INFTY) < EXP)
		{
			x.clear();
			z = INFTY;
			printf("unbounded");
			break;
		}
		PIVOT(max_thetai, min_ci,m,n);//做一次高斯消元
	}
	return pair<vector<double> ,double> (x, -z);
}
int INITLIZESIMPLEX(int m,int n)
{
	int min_bi = find_min_index(b_matrix, m);//找到bi最小的那个

	BI.clear();
	for (int i = m; i < n; i++)
	{
		BI.insert(i);
	}

	if (b_matrix[min_bi] >= 0)
		//左右的b都是非负的,那么这个问题就可以自然的获取一组顶点和解
	{
		return 1;
	}
	//使用黑科技,直接从bi最小的那一行 任选一个 小于0 的co_matrix[i]来操作
	int col_index = -1;
	for (int i = 0; i < n; i++)
	{
		if (co_matrix[min_bi][i] < 0)
		{
			col_index = i;
			break;
		}
	}
	if (col_index == -1)
	{
		BI.clear();
		printf("Infeasible");
		return 0;//无解
	}
	//否则执行一次转轴操作
	PIVOT(min_bi, col_index,m,n);
	return 1;
}
vector<double > CALCULATEX(int m,int n)
{
	vector<double> x;
	for (int i = 0; i < n; i++)
	{
		x.push_back(0);
	}
	for (auto it = BI.begin(); it != BI.end();it++)
	{
		int i = *it;
		{
			for (int j = 0; j < m; j++)
			{
				if (abs(co_matrix[j][i] -1 ) < EXP)
				{
					x[i] = b_matrix[j];
					if (x[i] < 0)
					{
						x.clear();
						printf("Infeasible");
						break;
						//无解
					}
				}
			}
		}
	}
	return x;
}
void PIVOT(int row, int col,int m,int n)
{
	//基的变换
	for (auto it = BI.begin(); it != BI.end(); it++)
		//从已有的基里面找到那个需要被替换出去的基
	{
		if (abs(co_matrix[row][*it] - 1) < EXP)
		{
			BI.erase(*it);
			break;
		}
	}
	BI.insert(col);

	//系数矩阵的高斯约旦消元
	for (int i = 0; i < n; i++)
		//把当前行除以co_matrix[row][col]
	{
		if (i != col)
		{
			co_matrix[row][i] /= co_matrix[row][col];
		}
	}
	b_matrix[row] /= co_matrix[row][col];
	co_matrix[row][col] = 1;

	//把其他行 减去第row行,校区第row行的第col个元素,以及b_matrix的高斯约旦
	for (int i = 0; i < m; i++)
	{
		if (i != row)
		{
			double factor = co_matrix[i][col] / co_matrix[row][col];
			for (int j = 0; j < n; j++)
			{
				co_matrix[i][j] -= co_matrix[row][j] * factor;
			}
			b_matrix[i] -= b_matrix[row] * factor;
		}
	}
	//c_matrix的高斯约旦
	double factor = c_matrix[col] / co_matrix[row][col];
	for (int j = 0; j < n; j++)
	{
		c_matrix[j] -= co_matrix[row][j] * factor;
	}
	//z的变换
	z -= b_matrix[row] * factor;



}
int main()
{
	//TEST 1.
	/*int m = 4, n = 7;
	co_matrix[0][0] = 1, co_matrix[0][1] = 1, co_matrix[0][2] = 1, co_matrix[0][3] = 1, co_matrix[0][4] = 0, co_matrix[0][5] = 0, co_matrix[0][6] = 0;
	co_matrix[1][0] = 1, co_matrix[1][1] = 0, co_matrix[1][2] = 0, co_matrix[1][3] = 0, co_matrix[1][4] = 1, co_matrix[1][5] = 0, co_matrix[1][6] = 0;
	co_matrix[2][0] = 0, co_matrix[2][1] = 0, co_matrix[2][2] = 1, co_matrix[2][3] = 0, co_matrix[2][4] = 0, co_matrix[2][5] = 1, co_matrix[2][6] = 0;
	co_matrix[3][0] = 0, co_matrix[3][1] = 3, co_matrix[3][2] = 1, co_matrix[3][3] = 0, co_matrix[3][4] = 0, co_matrix[3][5] = 0, co_matrix[3][6] = 1;
	b_matrix[0] = 4, b_matrix[1] = 2, b_matrix[2] = 3, b_matrix[3] = 6;
	c_matrix[0] = -1, c_matrix[1] = -14, c_matrix[2] = -6, c_matrix[3]= 0, c_matrix[4]= 0, c_matrix[5] = 0, c_matrix[6] = 0;*/

	//TEST 2.
	/*int m = 2, n = 4;
	co_matrix[0][0] = -1, co_matrix[0][1] = -1, co_matrix[0][2] = 1, co_matrix[0][3] = 0;
	co_matrix[1][0] = 1, co_matrix[1][1] = 1, co_matrix[1][2] = 0, co_matrix[1][3] = 1;
	b_matrix[0] = -1, b_matrix[1] = 2;
	c_matrix[0] = 1, c_matrix[1]= 2, c_matrix[2] = 0, c_matrix[3]= 0;
	*/
	//TEST 3.无解
	/*int m = 2, n = 4;
	co_matrix[0][0] = -1, co_matrix[0][1] = -1, co_matrix[0][2] = 1, co_matrix[0][3] = 0;
	co_matrix[1][0] = 1, co_matrix[1][1] = 1, co_matrix[1][2] = 0, co_matrix[1][3] = 1;
	b_matrix[0] = -2, b_matrix[1] = -1;
	c_matrix[0] = 1, c_matrix[1] = 2, c_matrix[2] = 0, c_matrix[3] = 0;*/
	//TEST 4.
	//int m = 3, n = 3;
	//co_matrix[0][0] = -1, co_matrix[0][1] = 0, co_matrix[0][2] = 0;
	//co_matrix[1][0] = -1, co_matrix[1][1] = -1, co_matrix[1][2] = 0 ;
	//co_matrix[2][0] = 0, co_matrix[2][1] = -1, co_matrix[2][2] = -1;
	//b_matrix[0] = -2, b_matrix[1] = -3,b_matrix[2]=-4;
	//c_matrix[0] = 2, c_matrix[1] = 5, c_matrix[2] = 2;
	int DAYN, VOLM;
		//天数	,志愿者的个数
	int m, n;
	scanf("%d%d", &DAYN, &VOLM);
	m = DAYN;//每天都是一个约束
	n = VOLM;//每种志愿者是一个变量
	for (int i = 0; i < DAYN; i++)
	//每天需要的志愿者个数
	{
		scanf("%f", &(b_matrix[i]));
		b_matrix[i] *= -1;
	}
	int si, fi;
	double ci;
	for (int i = 0; i <VOLM; i++)
	{
		scanf("%d%d%f", &si, &fi, &ci);
		for (int j = si - 1; j < fi; j++)
			//从si-1 到 fi行
		{
			co_matrix[j][i] = -1;
		}
		c_matrix[i] = ci;
	}
	auto rst = SIMPLEX(m, n);
	if (rst.first.size())
		printf("min:%0.5f\n", rst.second);

	system("pause");
	return 0;
}
