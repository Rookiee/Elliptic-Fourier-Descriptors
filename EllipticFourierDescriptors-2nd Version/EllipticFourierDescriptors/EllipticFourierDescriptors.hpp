/*
Created on May 29, 2016
@author: Haoyang Xie
C++ implementation of Elliptic Fourier Descriptors.
Original paper of the implemented algorithm is
"Elliptic Fourier of a Closed Contour", Frank P. Kuhl, Charles Giardian, 1981
Repository available at "https://github.com/Rookiee/Elliptic-Fourier-Descriptors"
Original MATLAB codes at "http://www.mathworks.com/matlabcentral/fileexchange/12746-elliptical-fourier-shape-descriptors/content/fEfourier.m"
and "http://www.mathworks.com/matlabcentral/fileexchange/12746-elliptical-fourier-shape-descriptors/content/rEfourier.m"

---Eigen library is necessary.
---"testdata.txt" in this directory is used to test this program.
---All of the data must be import to Eigen Matrix.
---If the data matrix is N x 2, use GetData2Cols function directly.
---If the data matrix is N x 3, use GetData3Cols firstly, then CreateMatrixWithoutYAxis function.
---The function fEfourier is used to compute the elliptic fourier descriptors of a contour.
---The function rEfourier is used to compute a reconstruct contour related to the input contour.

New features (2nd version):
	1. Design new input function "Eigen::MatrixXf getDataFromTxttoEigen(std::string fileName, int cols );" 
	   We now can get data from a .txt file directly by single function no matter 2 or 3 cols the file is.
	2. Add a output function "void saveDataFromEigentoTxt(std::string fileName, Eigen::MatrixXf& matrix);" 
	   This function could save a Eigen::MatrixXf into a .txt file.
    3. Add a new super function: "Eigen::MatrixXf getFinalFittingDataDirectly(std::string srcFileName, std::string dstFileName, int iNoOfHarmonics, bool bNormaliseSizeState, bool bNormaliseOrientationState,
																			int noOfPointsReconstruct, int cols, int deleteCol);"
	   The function could save the fitting file to a .txt file directly.what we need to do is just call it in our main() like:
	   			getFinalDataDirectly("testdata.txt", "testSave.txt", 10, true, true, 200, 3, 2);
	   				10 is the number of harmonics;
	   				200 is the points used to reconstruct the fitting shape
	   				3 is the number of cols in "testdata.txt"
	   				2 is the col removed if there are 3 cols in the "testdata.txt" 
	4. Separate the original single .h file into "EllipticFourierDescriptors.hpp" and "EllipticFourierDescriptors.cpp".

	
*/

#ifndef ELLIPTIC_FOURIER_DESCRIPTORS_H
#define ELLIPTIC_FOURIER_DESCRIPTORS_H


#include <vector>
#include <list>
#include <fstream>
#include <Eigen/Dense>
#define PI 3.1416

/* fileName: 
   cols: how many cols in the file. the default value is 3 */
Eigen::MatrixXf getDataFromTxttoEigen(std::string fileName, int cols );

/* if the file includes 3 cols, we need to delete a cols. 
   m: the matrix get from the above function.
   rows: the number of rows in the original Eigen::MatrixXf, the para could get from "tempMatrix.rows()" 
   cols: the number of cols is 3 usually.
   n: which col should be deleted. n is 1, 2 or 3 */
Eigen::MatrixXf CreateMatrixWithoutOneAxis(Eigen::MatrixXf m, int rows, int cols, int n);

/* 傅里叶级数展开:
Eigen::MatrixXf fEfourier(Eigen::MatrixXf outline, int iNoOfHarmonicsAnalyse, bool bNormaliseSizeState, bool bNormaliseOrientationState);
---MatrixXf outline(ROWS, 2),
---int iNoOfHarmonicsAnalyse: 谐波次数，
---bool bNormaliseSizeState
---bool bNormaliseOrientationState
****return: rFSDs-----MatrixXf(4,n)****/
Eigen::MatrixXf fEfourier(Eigen::MatrixXf outline, int iNoOfHarmonicsAnalyse, bool bNormaliseSizeState, bool bNormaliseOrientationState);

/* Get the fitting shape */
Eigen::MatrixXf rEfourier(Eigen::MatrixXf rFSDs, int iNoOfHarmonicsReconstruct, int iNoOfPointsReconstruct);

/* (Super Function) the function compute the final result directly. 
	cols: the number of cols in srcFileName. */
void getFinalFittingDataDirectly(std::string srcFileName, std::string dstFileName, int iNoOfHarmonics, bool bNormaliseSizeState, bool bNormaliseOrientationState,
	int noOfPointsReconstruct, int cols, int deleteCol);

void saveDataFromEigentoTxt(std::string fileName, Eigen::MatrixXf& matrix);

#endif