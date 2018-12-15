#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <vector>
#include <algorithm>
#include <set>
#include <utility>
#define EXP (0.2e-20)
using namespace std;
#define INFTY (2e20)
#pragma once
struct _Matrix
	//����Ķ���,�����ڱ�ʾϵ������,������,������
{
	double ** data;
	int cols, rows;
	char * name;
	_Matrix()
		:cols(0), rows(0), data(0), name("")
	{

	}
	void init()
	{
		data = (double **)malloc(sizeof(double*)* rows);
		for (int i = 0; i < rows; i++)
		{
			data[i] = (double*)malloc(sizeof(double)* cols);
			memset(data[i], 0, sizeof(double)* cols);
		}
	}
	_Matrix(int _rows, int _cols, char *_name = "Matrix")
		:cols(_cols), rows(_rows), name(_name), data(0)
	{
		init();
	}
	~_Matrix()
	{
		//if (data)
		//{
		//	for (int i = 0; i < rows; i++)
		//	{
		//		if (data[i])
		//			free(data[i]);
		//	}
		//	free(data);
		//}
		//rows = cols = 0;
	}
	_Matrix& operator= (_Matrix & rhs)
		//deep copy .���캯����ʹ��ǳ�㸴��
	{
		if (data != NULL)
		{
			this->~_Matrix();
		}
		rows = rhs.rows;
		cols = rhs.cols;
		name = rhs.name;
		init();
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				data[i][j] = rhs.data[i][j];
			}
		}
		return *this;
	}
	double* operator[](int rows)
	{
		return data[rows];
	}
	void display()
	{
		printf("Matrix of %s \n", name);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				printf("%0.5lf\t", data[i][j]);
			}
			printf("\n");
		}
	}
};
typedef _Matrix Mat;
typedef _Matrix * pMat;

class LP
{
public:
	LP()
	{

	}
	void set_cofficient_matrix(Mat & mat);
	//�����ɳ���
	void set_value_column(Mat & mat);
	void set_object_function(Mat & mat);

	pair< Mat, double >SIMPLEX();
	pair< set<int>, Mat>INIT_SOLUTION();//Ѱ��һ����ж���,�Լ���Ӧ��һ������е��±�
	void PIVOT(int i, int j);//ѡ��aijΪ��Ԫ,��aij���ڵĵ�j�����µ�Ԫȫ������

	~LP()
	{

	}
public:
	Mat co_matrix;
	Mat bi_matrix;
	Mat ci_matrix;
private:
	double z;//Ŀ�꺯��ֵ���෴��
	set<int> Bi;//һ������е��±�
	Mat X;//�����������¶���
};

