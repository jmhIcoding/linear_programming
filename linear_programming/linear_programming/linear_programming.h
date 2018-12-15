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
	//矩阵的定义,可用于表示系数矩阵,列向量,行向量
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
		//deep copy .构造函数就使用浅层复制
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
	//接受松弛型
	void set_value_column(Mat & mat);
	void set_object_function(Mat & mat);

	pair< Mat, double >SIMPLEX();
	pair< set<int>, Mat>INIT_SOLUTION();//寻找一组可行顶点,以及对应的一组基本列的下标
	void PIVOT(int i, int j);//选中aij为主元,把aij所在的第j列上下的元全部消掉

	~LP()
	{

	}
public:
	Mat co_matrix;
	Mat bi_matrix;
	Mat ci_matrix;
private:
	double z;//目标函数值的相反数
	set<int> Bi;//一组基本列的下标
	Mat X;//给定基本列下顶点
};

