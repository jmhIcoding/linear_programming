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
			//找到了最优解
		{
			//构建最优解
			return pair<Mat, double>(X, -z);
			break;
		}
		//说明还可以优化,那就选择这一列加入
		//计算非负的最小的theta
		Mat theta(this->co_matrix.rows, 1, "theta");
		//bi_matrix.display();

		for (int i = 0; i < this->co_matrix.rows; i++)
		{
			if (abs(co_matrix[i][index]) < EXP)
				//为0
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
			//不可解
			this->X.rows  = X.cols=0;
			this->X.init();
			this->z = 0;
			return pair<Mat, double>(X, z);
		}
		else if (abs(theta[rows_index][0] - INFTY)<EXP)
		{
			//没有最小值,没有下界
			return pair<Mat, double>(X, -INFTY);
		}
		//把(rows_index,index)作为主元,进行消元
		//这说明,接下来,将要把第index列换入,同时把

		
		//构造新的顶点
		//构造lambda
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
		//同时,把旧的基本列里面那些第_i行为1的列丢掉。
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
		//所有都是非负的,那么所有松弛变量都相应的bi
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
		//需要解辅助方程,待解决
	{

	}
}
void LP::PIVOT(int _i, int _j)//选中aij为主元,把aij所在的第j列上下的元全部消掉,同时自己消成1
{
	for (int j = 0; j < this->co_matrix.cols; j++)
		//把自己所在的行处理好,把主元化成1.
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
		//把主元所在列的其他元都消除
	{
		if (i != _i)
		{
			double factor = co_matrix[i][_j] / co_matrix[_i][_j];//求出系数,第i行的,第_j列与主元 第_i行,第_j列的熵
			for (int j = 0; j < this->co_matrix.cols; j++)
			{
				co_matrix[i][j] = co_matrix[i][j]  - co_matrix[_i][j] * factor;
			}
			//把bi处理一下
			bi_matrix[i][0] = bi_matrix[i][0] - bi_matrix[_i][0] * factor;
		}
	}
	//把ci也处理一下
	//co_matrix.display();
	//bi_matrix.display();
	double factor = ci_matrix[_j][0]/co_matrix[_i][_j];
	for (int i = 0; i < ci_matrix.rows; i++)
	{
		ci_matrix[i][0] -= factor * co_matrix[_i][i];
	}
	z -= factor * bi_matrix[_i][0];
	//ci_matrix.display();
	//主元是(_i,_j),将要把第_j 列换入.
	

}
