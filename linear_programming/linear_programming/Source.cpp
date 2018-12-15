#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <vector>
#include <algorithm>
#include "linear_programming.h"
using namespace std;
int main()
{
	//TEST 1.
	//Mat cofficient(4, 7, "cofficient");
	//Mat bi(4, 1, "bi");
	//Mat ci(7, 1, "ci");
	//cofficient[0][0] = 1, cofficient[0][1] = 1, cofficient[0][2] = 1, cofficient[0][3] = 1, cofficient[0][4] = 0, cofficient[0][5] = 0, cofficient[0][6] = 0;
	//cofficient[1][0] = 1, cofficient[1][1] = 0, cofficient[1][2] = 0, cofficient[1][3] = 0, cofficient[1][4] = 1, cofficient[1][5] = 0, cofficient[1][6] = 0;
	//cofficient[2][0] = 0, cofficient[2][1] = 0, cofficient[2][2] = 1, cofficient[2][3] = 0, cofficient[2][4] = 0, cofficient[2][5] = 1, cofficient[2][6] = 0;
	//cofficient[3][0] = 0, cofficient[3][1] = 3, cofficient[3][2] = 1, cofficient[3][3] = 0, cofficient[3][4] = 0, cofficient[3][5] = 0, cofficient[3][6] = 1;
	//bi[0][0] = 4, bi[1][0] = 2, bi[2][0] = 3, bi[3][0] = 6;
	//ci[0][0] = -1, ci[1][0] = -14, ci[2][0] = -6, ci[3][0] = 0, ci[4][0] = 0, ci[5][0] = 0, ci[6][0] = 0;
	//cofficient.display();
	//bi.display();
	//ci.display();
	//TEST 2.
	Mat cofficient(2, 4, "cofficient");
	Mat bi(2, 1, "bi");
	Mat ci(4, 1, "ci");
	cofficient[0][0] = -1, cofficient[0][1] = -1, cofficient[0][2] = 1, cofficient[0][3] = 0;
	cofficient[1][0] = 1, cofficient[1][1] = 1, cofficient[1][2] = 0, cofficient[1][3] = 1;
	bi[0][0] = -1, bi[1][0] = 2;
	ci[0][0] = 1, ci[1][0] = 2, ci[2][0] = 0, ci[3][0] = 0;
	LP test;
	test.set_cofficient_matrix(cofficient);
	test.set_object_function(ci);
	test.set_value_column(bi);
	auto rst=test.SIMPLEX();
	rst.first.display();
	printf("min:%0.5f\n", rst.second);
	system("pause");
	return 0;
}
