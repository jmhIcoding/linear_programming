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
								//aux ��־λ���ڿ���SIMPLEX�Ƿ������븨�����Թ滮����Ϊ�������Թ滮�� ��ʼ��Ƚ�����
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
		lambda.display();
		X.display();
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
	int l = 0;//��С��bi���±�
	for (int i = 0; i < bi_matrix.rows; i++)
	{
		if (min_bi> bi_matrix[i][0])
		{
			min_bi = bi_matrix[i][0];
			l = i;
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
		//���츨�����Թ滮
		Mat aux_co_matrix(co_matrix.rows, co_matrix.cols + 1, "aux cofficient matrix");
		Mat aux_bi_matrix(bi_matrix.rows, 1, "aux bi matrix");//m��Լ��
		Mat aux_ci_matrix(co_matrix.cols + 1, 1);//n+1��δ֪��
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
		//����ӻ�
		for (int i = aux_co_matrix.rows+1; i < aux_co_matrix.cols; i++)
		{
			if (aux_lp.co_matrix[l][i] != 1)
			{
				aux_lp.Bi.insert(i);
			}
		}
		//ͬʱ��x0���ڵ��е���һ��������
		aux_lp.Bi.insert(aux_co_matrix.rows);//�����û�
		//����x0���ڵ��еĵ�l��ԪΪ��Ԫ,��һ�θ�˹��Ԫ,
		aux_lp.PIVOT(l, aux_co_matrix.rows);
		//������������
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

		//�⸨�����Թ滮
		pair<Mat,double> aux_rst=aux_lp.SIMPLEX(1);
		if (abs(aux_rst.second) < EXP)
		//ԭ���Թ滮�н�
		{
			if (aux_lp.Bi.find(aux_lp.co_matrix.rows)!=aux_lp.Bi.end())
				aux_lp.Bi.erase(aux_lp.co_matrix.rows);//����еĻ���ɾ��x0��һ��,�϶���û�е�
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
			//�ٰ�ԭ���Թ滮��Լ���޸�
			this->bi_matrix = aux_lp.bi_matrix;
			//this->ci_matrix; ci�ȱ�����
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
			//����һ�θ�˹��Ԫ,��ԭ���Թ滮������ת
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
void LP::PIVOT(int _i, int _j)//ѡ��aijΪ��Ԫ,��aij���ڵĵ�j�����µ�Ԫȫ������,ͬʱ�Լ�����1
{
	for (int j = 0; j < this->co_matrix.cols; j++)
		//���Լ����ڵ��д����,����Ԫ����1.
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
	co_matrix.display();
	bi_matrix.display();
	double factor = ci_matrix[_j][0]/co_matrix[_i][_j];
	for (int i = 0; i < ci_matrix.rows; i++)
	{
		ci_matrix[i][0] -= factor * co_matrix[_i][i];
	}
	z -= factor * bi_matrix[_i][0];
	ci_matrix.display();
	//��Ԫ��(_i,_j),��Ҫ�ѵ�_j �л���.
	

}
