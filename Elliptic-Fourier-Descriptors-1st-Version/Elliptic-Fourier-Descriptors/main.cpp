
#include <QtCore/QCoreApplication>
#include "Elliptic_Fourier_Descriptors.h";

int main(int argc, char *argv[])
{
	QCoreApplication app(argc, argv);

	Eigen::MatrixXf m(ROWS, 2);
	m = GetData2Cols("testdata.txt");
	//std::cout << m << std::endl;
	/*Eigen::MatrixXf newM = CreateMatrixWithoutYAxis(m,2);*/
	//std::cout << newM << std::endl;
	
	const int noOfHarmonicsAnalyse = 10;
	Eigen::MatrixXf rFSDs(4, noOfHarmonicsAnalyse);
	rFSDs = fEfourier(m, noOfHarmonicsAnalyse, true, true);
	for (int i = 0; i < 4; ++i){
		rFSDs(i, 0) = 0;
	}

	const int noOfHarmonicsReconstruct = 10;
	const int noOfPointsReconstruct = 200;
	Eigen::MatrixXf outln(noOfPointsReconstruct, 2);
	outln = rEfourier(rFSDs, noOfHarmonicsReconstruct, noOfPointsReconstruct);
	std::cout << "      Reconstruct Result" <<std::endl;
	std::cout << "***********************************************" << std::endl;
	std::cout << outln << std::endl;
	std::cout << "            DONE" << std::endl;

	return app.exec();

	
}
