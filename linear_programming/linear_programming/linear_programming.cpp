#include "linear_programming.h"
void LP::set_cofficient_matrix(Mat & mat)
{
	this->co_matrix = mat;
	this->X.rows = mat.cols;
	this->X.cols = 1;
	this->X.init();
}
void LP::set_object_function(Mat & mat)
{
	this->ci_matrix = mat;
}
void LP::set_value_column(Mat & mat)
{
	this->bi_matrix = mat;
}
int find_min(Mat & mat)
{
	double min_ci = mat[0][0];
	int index = 0;
	for (int i = 0; i < mat.rows; i++)
	{
		if (min_ci> mat[i][0])
		{
			min_ci = mat[i][0];
			index = i;
		}
	}
	return index;
}
pair<Mat, double > LP::SIMPLEX()
{
	pair<set<int>, Mat > init = INIT_SOLUTION();
	this->Bi = init.first;
	this->X = init.second;
	while (true)
	{
		//find a columns to be added to the basic colums
		int index = find_min(this->ci_matrix);
		if (ci_matrix[index][0] >= 0)
			//�ҵ������Ž�
		{
			//�������Ž�
			return pair<Mat, double>(X, -z);
			break;
		}
		//˵���������Ż�,�Ǿ�ѡ����һ�м���
		//����Ǹ�����С��theta
		Mat theta(this->co_matrix.rows, 1, "theta");
		//bi_matrix.display();

		for (int i = 0; i < this->co_matrix.rows; i++)
		{
			if (abs(co_matrix[i][index]) < EXP)
				//Ϊ0
			{
				theta[i][0] = INFTY;
			}
			else
			{
				theta[i][0] = this->bi_matrix[i][0] / co_matrix[i][index];
			}
		}
		//theta.display();
		int rows_index = find_min(theta);
		if (theta[rows_index][0] < 0)
		{
			//���ɽ�
			this->X.rows  = X.cols=0;
			this->X.init();
			this->z = 0;
			return pair<Mat, double>(X, z);
		}
		else if (abs(theta[rows_index][0] - INFTY)<EXP)
		{
			//û����Сֵ,û���½�
			return pair<Mat, double>(X, -INFTY);
		}
		//��(rows_index,index)��Ϊ��Ԫ,������Ԫ
		//��˵��,������,��Ҫ�ѵ�index�л���,ͬʱ��

		
		//�����µĶ���
		//����lambda
		Mat lambda(this->co_matrix.cols, 1, "lambda");
		for (auto it = Bi.begin(); it != Bi.end(); it++)
		{
			for (int i = 0; i < co_matrix.rows; i++)
			{
				if (abs(co_matrix[i][*it] - 1) < EXP)
				{
					lambda[*it][0] = co_matrix[i][index];
				}
			}
		}
		lambda[index][0] = -1;//
		//lambda.display();
		//X.display();
		for (int i = 0; i < co_matrix.cols; i++)
		{
			X[i][0] = X[i][0] - theta[rows_index][0] * lambda[i][0];
		}
		//ͬʱ,�ѾɵĻ�����������Щ��_i��Ϊ1���ж�����
		//X.display();
		for (auto it = Bi.begin(); it != Bi.end(); it++)
		{
			if (abs(co_matrix[rows_index][*it] - 1) < EXP)
			{
				Bi.erase(it);
				break;
			}
		}
		Bi.insert(index);
		PIVOT(rows_index, index);
	}
}

pair< set<int>, Mat > LP::INIT_SOLUTION()
{
	set<int> basic_columns;
	Mat result(co_matrix.cols,1,"result");
	double min_bi = this->bi_matrix[0][0];
	for (int i = 0; i < bi_matrix.rows; i++)
	{
		if (min_bi> bi_matrix[i][0])
		{
			min_bi = bi_matrix[i][0];
		}
	}
	if (min_bi >= 0)
		//���ж��ǷǸ���,��ô�����ɳڱ�������Ӧ��bi
	{
		for (int i = 0; i < co_matrix.cols; i++)
		{
			if (i >=( this->co_matrix.rows-1))
			{

				basic_columns.insert(i);
				result[i][0]=this->bi_matrix[i- (this->co_matrix.rows-1)][0];
			}
			else
			{
				result[i][0]=0;
			}
		}
		return pair<set<int>,Mat>(basic_columns, result);
	}
	else
		//��Ҫ�⸨������,�����
	{

	}
}
void LP::PIVOT(int _i, int _j)//ѡ��aijΪ��Ԫ,��aij���ڵĵ�j�����µ�Ԫȫ������,ͬʱ�Լ�����1
{
	for (int j = 0; j < this->co_matrix.cols; j++)
		//���Լ����ڵ��д�����,����Ԫ����1.
	{
		if (j != _j)
		{
			co_matrix[_i][j] = co_matrix[_i][j] / co_matrix[_i][_j];
		}
	}
	bi_matrix[_i][0] /= co_matrix[_i][_j];
	co_matrix[_i][_j] = 1;
	//bi_matrix.display();
	//co_matrix.display();
	for (int i = 0; i < this->co_matrix.rows; i++)
		//����Ԫ�����е�����Ԫ������
	{
		if (i != _i)
		{
			double factor = co_matrix[i][_j] / co_matrix[_i][_j];//���ϵ��,��i�е�,��_j������Ԫ ��_i��,��_j�е���
			for (int j = 0; j < this->co_matrix.cols; j++)
			{
				co_matrix[i][j] = co_matrix[i][j]  - co_matrix[_i][j] * factor;
			}
			//��bi����һ��
			bi_matrix[i][0] = bi_matrix[i][0] - bi_matrix[_i][0] * factor;
		}
	}
	//��ciҲ����һ��
	//co_matrix.display();
	//bi_matrix.display();
	double factor = ci_matrix[_j][0]/co_matrix[_i][_j];
	for (int i = 0; i < ci_matrix.rows; i++)
	{
		ci_matrix[i][0] -= factor * co_matrix[_i][i];
	}
	z -= factor * bi_matrix[_i][0];
	//ci_matrix.display();
	//��Ԫ��(_i,_j),��Ҫ�ѵ�_j �л���.
	

}