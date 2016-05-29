/*
Created on May 24, 2016
@author: Haoyang Xie
C++ implementation of Elliptic Fourier Descriptors.
Original paper of the implemented algorithm is
"Elliptic Fourier of a Closed Contour", Frank P. Kuhl, Charles Giardian, 1981
Repository available at "https://github.com/Rookiee/Elliptic-Fourier-Descriptors"
Origianl MATLAB codes at "http://www.mathworks.com/matlabcentral/fileexchange/12746-elliptical-fourier-shape-descriptors/content/fEfourier.m"
and "http://www.mathworks.com/matlabcentral/fileexchange/12746-elliptical-fourier-shape-descriptors/content/rEfourier.m"

---Eigen library is necessary.
---"testdata.txt" in this directory is used to test this program.
---All of the data must be import to Eigen Matrix.
---If the data matrix is N x 2, use GetData2Cols function directly.
---If the data matrix is N x 3, use GetData3Cols firstly, then CreateMatrixWithoutYAxis function.
---The function fEfourier is used to compute the elliptic fourier descriptors of a contour.
---The function rEfourier is used to compute a reconstruct contour related to the input contour.
*/

#ifndef ELLIPTIC_FOURIER_DESCRIPTORS
#define ELLIPTIC_FOURIER_DESCRIPTORS


/* C Libraries */
#include <stdio.h>
#include <stdbool.h>
/* C++ Libraries */
#include <iostream>
#include <vector>
/* 3th Library : Eigen-3.2.8 */
#include "../include/eigen-eigen-3.2.8/Eigen/Dense"

/* 条件： 必须已知一个数据文件包含的行数和列数 */
#define ROWS 100
#define COLS 3
//#define COLS 2

#define pi 3.1416

/* 直接读取一个数据文件，并保存为Eigen::MatrixXf(ROWS, COLS)矩阵 */
Eigen::MatrixXf GetData3Cols(const char* filePath){
	Eigen::MatrixXf dataMatrix(ROWS, COLS);
	freopen(filePath, "r", stdin);
	for (int i = 0; i < ROWS; ++i){
		for (int j = 0; j < COLS; ++j){
			std::cin >> dataMatrix(i, j);
		}
	}
	return dataMatrix;
}
Eigen::MatrixXf GetData2Cols(const char* filePath){
	Eigen::MatrixXf dataMatrix(ROWS, 2);
	freopen(filePath, "r", stdin);
	for (int i = 0; i < ROWS; ++i){
		for (int j = 0; j < 2; ++j){
			std::cin >> dataMatrix(i, j);
		}
	}
	return dataMatrix;
}

/* 删除矩阵的第n列后组成组成新矩阵 */
Eigen::MatrixXf CreateMatrixWithoutYAxis(Eigen::MatrixXf m, int n = 1){
	if (n >= 3){
		std::cerr << "Fail." << std::endl;
		exit(EXIT_FAILURE);
	}
	Eigen::MatrixXf newM(ROWS, COLS - 1);
	if (n == 0){	//删除x,第一列

		newM.col(0) = m.col(COLS - 2);
		newM.col(1) = m.col(COLS - 1);
	}
	if (n == 1){	//删除y， 第二列
		newM.col(0) = m.col(0);
		newM.col(1) = m.col(COLS - 1);

	}
	if (n == 2){	//删除z， 第三列
		newM.col(0) = m.col(0);
		newM.col(1) = m.col(COLS - 2);

	}
	return newM;
}


/* 傅里叶级数展开:
Eigen::MatrixXf fEfourier(Eigen::MatrixXf outline, int iNoOfHarmonicsAnalyse, bool bNormaliseSizeState, bool bNormaliseOrientationState);
---MatrixXf outline(ROWS, 2),
---int iNoOfHarmonicsAnalyse: 谐波次数，
---bool bNormaliseSizeState
---bool bNormaliseOrientationState
****return: rFSDs-----MatrixXf(4,n)****
*/
Eigen::MatrixXf
fEfourier(Eigen::MatrixXf outline, int iNoOfHarmonicsAnalyse, bool bNormaliseSizeState, bool bNormaliseOrientationState){

	Eigen::MatrixXf rFSDs(4, iNoOfHarmonicsAnalyse);	// 最后返回的椭圆傅里叶算子

	/* 计算常量：2*pi*n */
	Eigen::MatrixXf iTwoNPi(1, iNoOfHarmonicsAnalyse);
	for (int i = 1; i <= iNoOfHarmonicsAnalyse; ++i){
		iTwoNPi(i - 1) = i;
	}
	iTwoNPi = iTwoNPi * 2 * pi;
	//std::cout << iTwoNPi << std::endl;


	/* 计算常亮：2*pi*pi*n*n */
	Eigen::MatrixXf rTwoNSqPiSq(1, iNoOfHarmonicsAnalyse);
	for (int i = 1; i <= iNoOfHarmonicsAnalyse; ++i){
		rTwoNSqPiSq(i - 1) = i*i;
	}
	rTwoNSqPiSq = rTwoNSqPiSq * 2 * pi*pi;
	//std::cout << rTwoNSqPiSq << std::endl;


	int iNoOfPoints = outline.rows() - 1; // hence there is 1 more data point in outline than iNoOfPoints

	Eigen::VectorXf rDeltaX = Eigen::VectorXf::Zero(iNoOfPoints + 1);	// X,Y,Z 列向量
	Eigen::VectorXf rDeltaY = Eigen::VectorXf::Zero(iNoOfPoints + 1);
	Eigen::VectorXf rDeltaT = Eigen::VectorXf::Zero(iNoOfPoints + 1);

	for (int iCount = 1; iCount <= iNoOfPoints; ++iCount){
		rDeltaX(iCount - 1) = outline(iCount, 0) - outline(iCount - 1, 0);
		rDeltaY(iCount - 1) = outline(iCount, 1) - outline(iCount - 1, 1);
	}
	//std::cout << rDeltaX << std::endl;
	//std::cout << rDeltaY << std::endl;

	for (int iCount = 0; iCount < iNoOfPoints + 1; ++iCount){
		rDeltaT(iCount) = sqrt((rDeltaX(iCount))*rDeltaX(iCount) + rDeltaY(iCount)*rDeltaY(iCount));
	}
	//std::cout << "-------------------" << std::endl << rDeltaT << std::endl;

	/* 删除X,Y,T中为0的元素 */
	std::vector<float> tempVector;
	int lenOfVector = 0;

	/* 对rDeltaT 操作*/
	for (int iCount = 0; iCount < iNoOfPoints; ++iCount){
		tempVector.push_back(rDeltaT(iCount));
	}
	lenOfVector = tempVector.size();
	Eigen::VectorXf tempEigenVector = Eigen::VectorXf::Zero(lenOfVector);
	for (int iCount = 0; iCount < lenOfVector; ++iCount){
		tempEigenVector(iCount) = tempVector[iCount];
	}
	rDeltaT = tempEigenVector;
	tempVector.clear();

	/*for (int iCount = 0; iCount < iNoOfPoints; ++iCount){
	if (rDeltaT(iCount) != 0){
	tempVector.push_back(rDeltaT(iCount));
	}
	}
	lenOfVector = tempVector.size();
	Eigen::VectorXf tempEigenVector = Eigen::VectorXf::Zero(lenOfVector);
	typedef std::vector<float>::iterator Iter;
	for (int iCount = 0; iCount < lenOfVector; ++iCount){
	tempEigenVector(iCount) = tempVector[iCount];
	}
	rDeltaT = tempEigenVector;
	tempVector.clear();*/

	//std::cout << "-------------------" << std::endl << rDeltaT << std::endl;

	/* 对deltaX操作*/
	for (int iCount = 0; iCount < iNoOfPoints; ++iCount){
		tempVector.push_back(rDeltaX(iCount));
	}
	lenOfVector = tempVector.size();
	tempEigenVector = Eigen::VectorXf::Zero(lenOfVector);
	for (int iCount = 0; iCount < lenOfVector; ++iCount){
		tempEigenVector(iCount) = tempVector[iCount];
	}
	rDeltaX = tempEigenVector;
	tempVector.clear();

	/*for (int iCount = 0; iCount < iNoOfPoints; ++iCount){
	if (rDeltaX(iCount) != 0){
	tempVector.push_back(rDeltaX(iCount));
	}
	}
	lenOfVector = tempVector.size();
	tempEigenVector = Eigen::VectorXf::Zero(lenOfVector);
	typedef std::vector<float>::iterator Iter;
	for (int iCount = 0; iCount < lenOfVector; ++iCount){
	tempEigenVector(iCount) = tempVector[iCount];
	}
	rDeltaX = tempEigenVector;
	tempVector.clear();*/
	//std::cout << "-------------------" << std::endl << rDeltaX << std::endl;
	/* 对deltaY操作*/
	/* 根据Matlab代码，此操作只是删除向量末尾的0，对于中间的0不做处理
	Matlab代码注释不清楚。*/
	for (int iCount = 0; iCount < iNoOfPoints; ++iCount){
		tempVector.push_back(rDeltaY(iCount));
	}
	lenOfVector = tempVector.size();
	tempEigenVector = Eigen::VectorXf::Zero(lenOfVector);
	for (int iCount = 0; iCount < lenOfVector; ++iCount){
		tempEigenVector(iCount) = tempVector[iCount];
	}
	rDeltaY = tempEigenVector;
	tempVector.clear();
	//std::cout << "BBBBBBBBBBBBBBBBBB" << std::endl;
	//std::cout << rDeltaY.size() << std::endl;
	//std::cout << rDeltaY << std::endl;

	/*
	for (int iCount = 0; iCount < iNoOfPoints; ++iCount){
	if (rDeltaY(iCount) != 0){
	tempVector.push_back(rDeltaY(iCount));
	}
	}
	lenOfVector = tempVector.size();
	tempEigenVector = Eigen::VectorXf::Zero(lenOfVector);
	typedef std::vector<float>::iterator Iter;
	for (int iCount = 0; iCount < lenOfVector; ++iCount){
	tempEigenVector(iCount) = tempVector[iCount];
	}
	rDeltaY = tempEigenVector;
	tempVector.clear();*/



	/* 减去一个重复点*/
	iNoOfPoints = rDeltaT.size() - 1;


	Eigen::VectorXf rTime = Eigen::VectorXf::Zero(iNoOfPoints + 1);

	for (int iCount = 1; iCount <= iNoOfPoints; ++iCount){
		rTime(iCount) = rTime(iCount - 1) + rDeltaT(iCount - 1);
	}
	//std::cout << "______________________________________________________" << std::endl;
	//std::cout << rTime.transpose() << std::endl;
	float rPeriod = rTime(iNoOfPoints);

	float rSum1 = 0;
	float rInnerDiff = 0;
	/* calculate the A-sub-0 coefficient*/
	for (int iP = 1; iP < iNoOfPoints + 1; ++iP){
		float rSum2 = 0;
		float rSum3 = 0;
		float rInnerDiff = 0;
		if (iP > 0){
			for (int iJ = 1; iJ < iP; ++iJ){
				rSum2 = rSum2 + rDeltaX(iJ - 1);
				rSum3 = rSum3 + rDeltaT(iJ - 1);
			}
			rInnerDiff = rSum2 - ((rDeltaX(iP - 1) / rDeltaT(iP - 1))*rSum3);
		}
		float rIncr1 = ((rDeltaX(iP - 1) / (2 * rDeltaT(iP - 1)))* ((rTime(iP)*rTime(iP)) - rTime(iP - 1)*rTime(iP - 1)) + rInnerDiff*(rTime(iP) - rTime(iP - 1)));
		rSum1 += rIncr1;
	}
	//std::cout << "(((((((((((((((" << std::endl;
	//std::cout << rSum1 << std::endl;



	rFSDs(0, 0) = ((1 / rPeriod) * rSum1) + outline(0, 0);
	//std::cout << rFSDs(0, 0) << std::endl;
	/* 计算an */
	for (int iHNo = 1; iHNo < iNoOfHarmonicsAnalyse; ++iHNo){
		rSum1 = 0;
		for (int iP = 0; iP < iNoOfPoints; ++iP){
			float rIncr1 = (rDeltaX(iP) / rDeltaT(iP)) * ((cos(iTwoNPi(iHNo - 1)*rTime(iP + 1) / rPeriod) - cos(iTwoNPi(iHNo - 1)*rTime(iP) / rPeriod)));
			rSum1 += rIncr1;
		}
		rFSDs(0, iHNo) = (rPeriod / rTwoNSqPiSq(iHNo - 1)) * rSum1;
	}
	//std::cout << "+++++++++++++++++" << rSum1 << std::endl;
	//std::cout << rFSDs.row(0) << std::endl;


	rFSDs(1, 0) = 0;
	/* 计算bn */
	for (int iHNo = 1; iHNo < iNoOfHarmonicsAnalyse; ++iHNo){
		rSum1 = 0;
		for (int iP = 0; iP < iNoOfPoints; ++iP){
			float rIncr1 = (rDeltaX(iP) / rDeltaT(iP)) * ((sin(iTwoNPi(iHNo - 1)*rTime(iP + 1) / rPeriod) - sin(iTwoNPi(iHNo - 1)*rTime(iP) / rPeriod)));
			rSum1 += rIncr1;
		}
		rFSDs(1, iHNo) = (rPeriod / rTwoNSqPiSq(iHNo - 1)) * rSum1;
	}
	//std::cout << "###################" << std::endl;
	//std::cout << rFSDs.row(1) << std::endl;

	/* Calculate the C-sub-0 coefficient */
	rSum1 = 0;
	for (int iP = 1; iP < iNoOfPoints + 1; ++iP){
		float rSum2 = 0;
		float rSum3 = 0;
		float rInnerDiff = 0;
		/* calculate the partial sums - there are 0 for iCount = 1 */
		if (iP > 0){
			for (int iJ = 1; iJ < iP; ++iJ){
				rSum2 += rDeltaY(iJ - 1);
				rSum3 += rDeltaT(iJ - 1);

			}
			rInnerDiff = rSum2 - ((rDeltaY(iP - 1) / rDeltaT(iP - 1)) * rSum3);
		}
		float rIncr1 = ((rDeltaY(iP - 1) / (2 * rDeltaT(iP - 1)))* ((rTime(iP)*rTime(iP)) - rTime(iP - 1)*rTime(iP - 1)) + rInnerDiff*(rTime(iP) - rTime(iP - 1)));
		rSum1 += rIncr1;
	}
	rFSDs(2, 0) = ((1 / rPeriod)*rSum1) + outline(0, 1);
	//std::cout << "CCCCCCCCCCCCCCCC" << std::endl;
	//std::cout << rFSDs.row(2) << std::endl;



	/* calculate the C-sub-n */
	for (int iHNo = 1; iHNo < iNoOfHarmonicsAnalyse; ++iHNo){
		rSum1 = 0;
		for (int iP = 0; iP < iNoOfPoints; ++iP){
			float rIncr1 = (rDeltaY(iP) / rDeltaT(iP)) * ((cos(iTwoNPi(iHNo - 1)*rTime(iP + 1) / rPeriod) - cos(iTwoNPi(iHNo - 1)*rTime(iP) / rPeriod)));
			rSum1 += rIncr1;
		}
		rFSDs(2, iHNo) = (rPeriod / rTwoNSqPiSq(iHNo - 1)) * rSum1;
	}
	//std::cout << "DDDDDDDDDDDDDDD" << std::endl;
	//std::cout << rFSDs.row(2) << std::endl;
	rFSDs(3, 0) = 0;
	/* calculate the D-sub-n */
	for (int iHNo = 1; iHNo < iNoOfHarmonicsAnalyse; ++iHNo){
		rSum1 = 0;
		for (int iP = 0; iP < iNoOfPoints; ++iP){
			float rIncr1 = (rDeltaY(iP) / rDeltaT(iP)) * ((sin(iTwoNPi(iHNo - 1)*rTime(iP + 1) / rPeriod) - sin(iTwoNPi(iHNo - 1)*rTime(iP) / rPeriod)));
			rSum1 += rIncr1;
		}
		rFSDs(3, iHNo) = (rPeriod / rTwoNSqPiSq(iHNo - 1)) * rSum1;
	}
	//std::cout << "EEEEEEEEEEEEEEEEE" << std::endl;
	//std::cout << rFSDs.row(3) << std::endl;



	/* 标准化 */
	if (bNormaliseSizeState == true || bNormaliseOrientationState == true){
		Eigen::MatrixXf rFSDsTemp(rFSDs);
		float rTheta1 = 0.5 * atan(2 * (rFSDsTemp(0, 1) * rFSDsTemp(1, 1) + rFSDsTemp(2, 1) * rFSDsTemp(3, 1)) /
			(rFSDsTemp(0, 1)*rFSDsTemp(0, 1) + rFSDsTemp(2, 1)*rFSDsTemp(2, 1) - rFSDsTemp(1, 1)*rFSDsTemp(1, 1) - rFSDsTemp(3, 1)*rFSDsTemp(3, 1)));
		Eigen::MatrixXf rStarFSDs(4, iNoOfHarmonicsAnalyse);
		for (int iHNo = 0; iHNo < iNoOfHarmonicsAnalyse; ++iHNo){
			rStarFSDs(0, iHNo) = cos((iHNo)*rTheta1) * rFSDsTemp(0, iHNo) + sin(iHNo*rTheta1)*rFSDsTemp(1, iHNo);
			rStarFSDs(1, iHNo) = -sin((iHNo)*rTheta1) * rFSDsTemp(0, iHNo) + cos(iHNo*rTheta1)*rFSDsTemp(1, iHNo);
			rStarFSDs(2, iHNo) = cos((iHNo)* rTheta1) * rFSDsTemp(2, iHNo) + sin((iHNo)* rTheta1) * rFSDsTemp(3, iHNo);
			rStarFSDs(3, iHNo) = -sin((iHNo)* rTheta1) * rFSDsTemp(2, iHNo) + cos((iHNo)* rTheta1) * rFSDsTemp(3, iHNo);

		}

		float rPsi1 = atan(rStarFSDs(2, 1) / rStarFSDs(0, 1));
		float rSemiMajor = sqrt(rStarFSDs(0, 1)*rStarFSDs(0, 1) + rStarFSDs(2, 1) * rStarFSDs(2, 1));
		rFSDs = rStarFSDs / rSemiMajor;

		//std::cout << "NNNNNNNNNNNNNN" << std::endl <<  rFSDs << std::endl;

		if (bNormaliseOrientationState == true){
			for (int iHNo = 0; iHNo < iNoOfHarmonicsAnalyse; ++iHNo){
				rFSDsTemp(0, iHNo) = (cos(rPsi1) * rStarFSDs(0, iHNo) + sin(rPsi1)* rStarFSDs(2, iHNo)) / rSemiMajor;
				rFSDsTemp(1, iHNo) = (cos(rPsi1) * rStarFSDs(1, iHNo) + sin(rPsi1)* rStarFSDs(3, iHNo)) / rSemiMajor;
				rFSDsTemp(2, iHNo) = (-sin(rPsi1) * rStarFSDs(0, iHNo) + cos(rPsi1)* rStarFSDs(2, iHNo)) / rSemiMajor;
				rFSDsTemp(3, iHNo) = (-sin(rPsi1) * rStarFSDs(1, iHNo) + cos(rPsi1)*rStarFSDs(3, iHNo)) / rSemiMajor;

			}
			rFSDs = rFSDsTemp;
			//std::cout << "MMMMMMMMMMMMMMM\n" << rFSDs << std::endl;
		}
	}

	return rFSDs;
}

Eigen::MatrixXf
rEfourier(Eigen::MatrixXf rFSDs, int iNoOfHarmonicsReconstruct, int iNoOfPointsReconstruct){
	Eigen::MatrixXf outln(iNoOfPointsReconstruct, 2);
	int iStartHarmonic = 1;
	Eigen::MatrixXf ReconnedOutline(iNoOfPointsReconstruct, 2);

	//std::cout << "***************\n " << rFSDs << "\n" << "********************" << std::endl;
	/* reconstruct the x-projection */
	for (int iTime = 0; iTime < iNoOfPointsReconstruct; ++iTime){
		float rSum = 0.0;
		for (int iHNo = iStartHarmonic; iHNo < iNoOfHarmonicsReconstruct; ++iHNo){
			rSum += (rFSDs(0, iHNo) * cos(2 * (iHNo)*pi*(iTime + 1) / iNoOfPointsReconstruct) +
				rFSDs(1, iHNo) * sin(2 * (iHNo)*pi*(iTime + 1) / iNoOfPointsReconstruct));
		}
		ReconnedOutline(iTime, 0) = rFSDs(0, 0) + rSum;
	}

	//std::cout << "QQQQQQQQQQQQQQ\n" << ReconnedOutline << std::endl;

	/* reconstruct the y-projection */
	for (int iTime = 0; iTime < iNoOfPointsReconstruct; ++iTime){
		float rSum = 0.0;
		for (int iHNo = iStartHarmonic; iHNo < iNoOfHarmonicsReconstruct; ++iHNo){
			rSum += (rFSDs(2, iHNo) * cos(2 * (iHNo)*pi*(iTime + 1) / iNoOfPointsReconstruct) +
				rFSDs(3, iHNo)* sin(2 * (iHNo)*pi*(iTime + 1) / iNoOfPointsReconstruct));
		}
		ReconnedOutline(iTime, 1) = rFSDs(2, 0) + rSum;
	}
	outln = ReconnedOutline;
	//std::cout << "***************\n " << outln << "\n" << "********************" << std::endl;
	return outln;
}

#endif