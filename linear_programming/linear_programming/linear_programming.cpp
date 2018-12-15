#include "linear_programming.h"
void LP::set_cofficient_matrix(Mat & mat)
{
	this->co_matrix = mat;
	this->X.rows = mat.cols;
	this->X.cols = 1;
	this->X.init();
	this->z = 0;
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
pair<Mat, double > LP::SIMPLEX(int aux)
								//aux 标志位用于控制SIMPLEX是否是用与辅助线性规划。因为辅助线性规划的 初始解比较特殊
{
	if (aux == 0)
	{
		pair<set<int>, Mat > init = INIT_SOLUTION();
		this->Bi = init.first;
		this->X = init.second;
	}
	while (true)
	{
		//find a columns to be added to the basic colums
		ci_matrix.display();
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
		lambda.display();
		X.display();
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
	int l = 0;//最小的bi的下标
	for (int i = 0; i < bi_matrix.rows; i++)
	{
		if (min_bi> bi_matrix[i][0])
		{
			min_bi = bi_matrix[i][0];
			l = i;
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
		//构造辅助线性规划
		Mat aux_co_matrix(co_matrix.rows, co_matrix.cols + 1, "aux cofficient matrix");
		Mat aux_bi_matrix(bi_matrix.rows, 1, "aux bi matrix");//m个约束
		Mat aux_ci_matrix(co_matrix.cols + 1, 1);//n+1个未知数
		aux_bi_matrix = bi_matrix;//deep copy.
		aux_ci_matrix[co_matrix.rows][0] = 1;
		for (int i = 0; i < aux_co_matrix.rows; i++)
		{
			for (int j = 0; j < aux_co_matrix.cols; j++)
			{
				if (j < co_matrix.rows)
				{
					aux_co_matrix[i][j] = co_matrix[i][j];
				}
				else if (j == co_matrix.rows)
				{
					aux_co_matrix[i][j] = -1;
				}
				else
				{
					aux_co_matrix[i][j] = co_matrix[i][j - 1];
				}
			}
		}
		LP aux_lp;
		aux_lp.set_cofficient_matrix(aux_co_matrix);
		aux_lp.set_object_function(aux_ci_matrix);
		aux_lp.set_value_column(aux_bi_matrix);
		//aux_lp.co_matrix.display();
		//aux_lp.ci_matrix.display();
		//aux_lp.bi_matrix.display();
		//先添加基
		for (int i = aux_co_matrix.rows+1; i < aux_co_matrix.cols; i++)
		{
			if (aux_lp.co_matrix[l][i] != 1)
			{
				aux_lp.Bi.insert(i);
			}
		}
		//同时把x0所在的列当做一个基本列
		aux_lp.Bi.insert(aux_co_matrix.rows);//构建好基
		//先以x0所在的列的第l个元为主元,做一次高斯消元,
		aux_lp.PIVOT(l, aux_co_matrix.rows);
		//接下来构建解
		aux_lp.co_matrix.display();
		aux_lp.ci_matrix.display();
		aux_lp.bi_matrix.display();
		for (int i = 0; i < aux_lp.co_matrix.cols; i++)
		{
			if (aux_lp.Bi.find(i) == aux_lp.Bi.end())
			{

				aux_lp.X[i][0] = 0;
			}
			else
			{
				int j = 0;
				for (j = 0; j < aux_lp.co_matrix.rows; j++)
				{
					if (abs(aux_lp.co_matrix[j][i] - 1) < EXP)
					{
						break;
					}
				}
				aux_lp.X[i][0]=aux_lp.bi_matrix[j][0];
			}
		}
		aux_lp.X.display();
		aux_lp.co_matrix.display();
		aux_lp.ci_matrix.display();
		aux_lp.bi_matrix.display();

		//解辅助线性规划
		pair<Mat,double> aux_rst=aux_lp.SIMPLEX(1);
		if (abs(aux_rst.second) < EXP)
		//原线性规划有解
		{
			if (aux_lp.Bi.find(aux_lp.co_matrix.rows)!=aux_lp.Bi.end())
				aux_lp.Bi.erase(aux_lp.co_matrix.rows);//如果有的话就删除x0那一列,肯定是没有的
			Bi = aux_lp.Bi;
			aux_lp.X.display();
			for (int i = 0; i < aux_lp.co_matrix.cols; i++)
			{
				if (i < this->co_matrix.rows)
				{

					this->X[i][0] = aux_lp.X[i][0];
				}
				else if (i>this->co_matrix.rows)
				{
					this->X[i-1][0] = aux_lp.X[i][0];
				}
			}
			//再把原线性规划的约束修改
			this->bi_matrix = aux_lp.bi_matrix;
			//this->ci_matrix; ci先保持着
			for (int i = 0; i < co_matrix.rows; i++)
			{
				for (int j = 0; j < aux_lp.co_matrix.cols; j++)
				{
					if (j < co_matrix.rows)
					{
						this->co_matrix[i][j] = aux_lp.co_matrix[i][j];
					}
					else if (j>co_matrix.rows)
					{
						this->co_matrix[i][j-1] = aux_lp.co_matrix[i][j];
					}
				}
			}
			//再做一次高斯消元,让原线性规划得以运转
			//int i, j;
			//for (auto it = Bi.begin(); it != Bi.end(); it++)
			//{
			//	if (*it < co_matrix.rows && abs(ci_matrix[*it][0])>EXP)
			//	{
			//		j = *it;
			//	}
			//}
			//for ( i = 0; i < co_matrix.rows; i++)
			//{
			//	if (abs(co_matrix[i][j])>EXP)
			//	{
			//		break;
			//	}
			//}
			//PIVOT(i, j);
			for (int i = 0; i < ci_matrix.rows; i++)
			{
				ci_matrix[i][0] *= -1;
			}
			return pair<set<int>, Mat>(aux_lp.Bi, X);
		}
		else
		{
			printf("unconsistent \n");
			Bi.clear();
			X.cols = 0;
			X.rows = 0;
			X.init();
			return pair< set<int>, Mat >(Bi,X);
		}

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
	co_matrix.display();
	bi_matrix.display();
	double factor = ci_matrix[_j][0]/co_matrix[_i][_j];
	for (int i = 0; i < ci_matrix.rows; i++)
	{
		ci_matrix[i][0] -= factor * co_matrix[_i][i];
	}
	z -= factor * bi_matrix[_i][0];
	ci_matrix.display();
	//主元是(_i,_j),将要把第_j 列换入.
	

}
