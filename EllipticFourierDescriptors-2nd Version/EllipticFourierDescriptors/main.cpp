#include <iostream>
#include "EllipticFourierDescriptors.hpp"

using namespace std;
int main(void){

	//Eigen::MatrixXf matrix;
	//// get data from a file with 3 cols.
	//matrix = getDataFromTxttoEigen("testdata.txt", 3);
	//// if there are three cols, we need to get the rows which will be used as 
	//// a parameter in next fuction to delete a 
	//int rows = matrix.rows();

	//matrix = CreateMatrixWithoutOneAxis(matrix, rows, 3, 3);
	//

	//Eigen::MatrixXf FDs = fEfourier(matrix, 10, true, true);
	//
	//for (int i = 0; i < 4; ++i){
	//	FDs(i, 0) = 0;
	//}

	//Eigen::MatrixXf outln = rEfourier(FDs, 10, 200);
	//
	//saveDataFromEigentoTxt("saveData.txt", outln);
	//****************************************************************************
	getFinalFittingDataDirectly("testdata.txt", "testSave.txt", 10, true, true, 100, 3, 2);
	system("pause");
	return 0;
}