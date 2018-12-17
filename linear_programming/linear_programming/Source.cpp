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
//INFIY ������
//EXP	����С��,��ҪΪ�˶Ը�double ���ж�ʱ���ܳ������
double co_matrix[N][N] = { 0 };
double b_matrix[N]= { 0 };
double c_matrix[N] = { 0 };
set<int> BI;//�洢�����е��±�
double z;//�洢��Сֵ
double theta[N] = { 0 };
double lambda[N] = { 0 };
int INITLIZESIMPLEX(int m,int n);
vector<double > CALCULATEX(int m,int n);
void PIVOT(int i,int j,int m,int n);//��co_matrix��aijΪ��Ԫ,ִ��һ�θ�˹Լ����Ԫ
int find_min_index(double * vec,int len)
{
	double tmp = vec[0];
	int index = 0;
	for (int i = 0; i < len; i++)
	{
		if (tmp>vec[i])
		{
			tmp = vec[i], index = i;//ѡ��1����С��
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
			tmp = vec[i], index = i;//ѡ��1������
		}
	}
	return index;
}
pair < vector<double>, double > SIMPLEX(int m, int n, int init_flag = 0)
{
	vector<double> x;
	int feasible = INITLIZESIMPLEX(m, n);//�����н��ʱ��
	while (feasible)
	{
		int min_ci = find_min_index(c_matrix, n);
		if (c_matrix[min_ci] >= 0)
			//���е�ci������0,˵��Ŀ�꺯���Ѿ������� �ҵ�������ֵ
		{
			x = CALCULATEX(m,n);
			break;
		}
		//��Ե�min_ci�� ������С�ķǸ�theta
		for (int i = 0; i < m; i++)
		{
			if (co_matrix[i][min_ci]>0)
				//���������theta���ǷǸ���
			{
				theta[i] = b_matrix[i] / co_matrix[i][min_ci];
			}
			else
				//��co_matrix[i][min_ci]Ϊ0,����Ϊ����ʱ��,������ȥ����
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
		PIVOT(max_thetai, min_ci,m,n);//��һ�θ�˹��Ԫ
	}
	return pair<vector<double> ,double> (x, -z);
}
int INITLIZESIMPLEX(int m,int n)
{
	int min_bi = find_min_index(b_matrix, m);//�ҵ�bi��С���Ǹ�

	BI.clear();
	for (int i = m; i < n; i++)
	{
		BI.insert(i);
	}

	if (b_matrix[min_bi] >= 0)
		//���ҵ�b���ǷǸ���,��ô�������Ϳ�����Ȼ�Ļ�ȡһ�鶥��ͽ�
	{
		return 1;
	}
	//ʹ�úڿƼ�,ֱ�Ӵ�bi��С����һ�� ��ѡһ�� С��0 ��co_matrix[i]������
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
		return 0;//�޽�
	}
	//����ִ��һ��ת�����
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
						//�޽�
					}
				}
			}
		}
	}
	return x;
}
void PIVOT(int row, int col,int m,int n)
{
	//���ı任
	for (auto it = BI.begin(); it != BI.end(); it++)
		//�����еĻ������ҵ��Ǹ���Ҫ���滻��ȥ�Ļ�
	{
		if (abs(co_matrix[row][*it] - 1) < EXP)
		{
			BI.erase(*it);
			break;
		}
	}
	BI.insert(col);

	//ϵ������ĸ�˹Լ����Ԫ
	for (int i = 0; i < n; i++)
		//�ѵ�ǰ�г���co_matrix[row][col]
	{
		if (i != col)
		{
			co_matrix[row][i] /= co_matrix[row][col];
		}
	}
	b_matrix[row] /= co_matrix[row][col];
	co_matrix[row][col] = 1;

	//�������� ��ȥ��row��,У����row�еĵ�col��Ԫ��,�Լ�b_matrix�ĸ�˹Լ��
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
	//c_matrix�ĸ�˹Լ��
	double factor = c_matrix[col] / co_matrix[row][col];
	for (int j = 0; j < n; j++)
	{
		c_matrix[j] -= co_matrix[row][j] * factor;
	}
	//z�ı任
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
	//TEST 3.�޽�
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
		//����	,־Ը�ߵĸ���
	int m, n;
	scanf("%d%d", &DAYN, &VOLM);
	m = DAYN;//ÿ�춼��һ��Լ��
	n = VOLM;//ÿ��־Ը����һ������
	for (int i = 0; i < DAYN; i++)
	//ÿ����Ҫ��־Ը�߸���
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
			//��si-1 �� fi��
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
