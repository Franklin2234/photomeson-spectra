//
// Created by Sergej on 20.06.18.
//
#include "calculator.h"
#include "ResultList.h"

/* ==================================================================================================================================================================== */
/* ======================================================================  GLOBAL VARIABLES  ========================================================================== */
/* ==================================================================================================================================================================== */


/* ================= Photon distribution ================ */

double alphaI = -1.6;
double alphaII = -1.8;
double N0 = 1;
double photonEnergy = 0.0;                          // [eV]






													/* ================= Proton distribution ================ */

double s = -2;
double K0 = 1;
double protonEnergy = 0.0;






/* ================ interaction variables =============== */

double CHI_pi0LR = 0.22;
double CHI_pi0HR = 0.26;
double CHI_piPlusLR = 0.22;
double CHI_piPlusHR = 0.25;
double CHI_Nu = 5.1E-4;

double M_LR_pi0 = 2.0 / 3.0;
double M_HR_pi0 = 0.47;
double M_piPlusLR = 1.0 / 3.0;
double M_piPlusHR = 0.77;

double sigmaLR = 2.0E-32;
double sigmaHR = 0.9E-32;







/* ========== repetitioneric variables ========= */

long int count = 0;
double y = 0.0;
double pionEnergy = 0.0;
double valueQ_pi = 0.0;
double valueResponse = 0.0;
double valuePhotonn = 0.0;
double valueProtonn = 0.0;
double lowBound = 0.5*145E6;
double upBound = 1E19;


std::vector<double> protonGrid;
std::vector<double> Q_pi;



/* ==================================================================================================================================================================== */
/* ======================================================================  calculatorULATIONS  ============================================================================== */
/* ==================================================================================================================================================================== */


double calculator::photonDistPrep(double y)
{
	photonEnergy = (M_PROTON*y) / protonEnergy;

	if (1E-3 <= photonEnergy && photonEnergy < 140)
	{
		valuePhotonn = N0 * pow((M_PROTON * CHI_pi0LR * photonEnergy) / (1E9*pionEnergy), alphaI);
	}
	else if (140 <= photonEnergy && photonEnergy < 3.6E3)
	{
		valuePhotonn = N0 * pow(1.4E-7, 0.2) * pow((M_PROTON * CHI_pi0LR * photonEnergy) / (1E9*pionEnergy), alphaII);
	}
	else
	{
		return 0.0;
	}
	return valuePhotonn;
}


double calculator::photonDistI(double x)
{
	valuePhotonn = N0 * pow((M_PROTON * CHI_pi0LR * x) / 1E9*pionEnergy, alphaI);
	return valuePhotonn;
}


double calculator::photonDistII(double x)
{
	valuePhotonn = N0 * pow((M_PROTON * CHI_pi0LR * x) / 1E9*pionEnergy, alphaII);
	return valuePhotonn;
}


double calculator::protonDistPrep(double x)
{
	if (1E9 <= x)
	{
		valueProtonn = K0 * pow(x / 1E9, s) * exp(-pow(x / 6.9E17, 2));
	}
	else
	{
		valueProtonn = 0.0;
	}
	return valueProtonn;
}


// Es wird für eine Photonenenergie die y-Werte bestimmt
void calculator::yList()
{
	std::ofstream listY;
	listY.open("Ylist_25.07.txt");

	std::ifstream readFile("listProtonEnergy.txt");
	while (readFile >> protonEnergy)
	{
		protonGrid.push_back(protonEnergy);             // die Protonenenergie wurde manuell in eine Datei geschrieben die hier eingelesen wird
	}                                                   // um immer das gleiche Grid der Pionenenergie zu bekommen. Danach wird es in die Ylist.txt geschrieben und in calculatorulate wieder aufgerufen

	for (photonEnergy = 1E-2; photonEnergy <= 3.6E3; photonEnergy += 1)
	{
		for (int i = 0; i < protonGrid.size(); i++)
		{
			y = (protonGrid[i] * photonEnergy) / M_PROTON;
			pionEnergy = CHI_pi0LR * protonGrid[i];
			listY << y << " " << photonEnergy << " " << protonGrid[i] << " " << pionEnergy << std::endl;
			std::cout << "Inner loop is runnung " << 100 * protonGrid[i] / upBound << " %." << std::endl;
			count++;
		}
		std::cout << "Outer loop is running. (y-Values) " << (100 * photonEnergy) / 3.6E3 << " %." << std::endl;
	}
	listY.close();
}


double calculator::responseFunc(double y) // Hümmer 31 und 32
{
	if (2 * y < 0.2E9)
	{
		valueResponse = 0.0;
	}
	if (0.2E9 <= 2 * y && y < 0.5E9)
	{
		valueResponse = M_LR_pi0 * sigmaLR*(1 - pow((0.2E9) / (2.0*y*y), 2.0));
	}
	if (0.5E9 <= 2 * y && y < 1.2E9)
	{
		valueResponse = M_HR_pi0 * sigmaHR*(1 - pow((0.5E9) / (2.0*y*y), 2.0));
	}
	if (1.2E9 <= 2 * y)
	{
		valueResponse = M_HR_pi0 * sigmaHR*((pow(1.2E9, 2) - pow(0.5E9, 2.0)) / pow(2 * y*y, 2.0));
	}
	return valueResponse;
}


double calculator::responseFunc2(double y) // Hümmer 31 und 32
{
	if (2 * y < 0.2E9) { valueResponse = 0.0; }
	if (0.2E9 <= 2 * y && y < 0.5E9) { valueResponse = M_LR_pi0 * sigmaLR*(1 - pow((0.2E9) / (2.0*y), 2.0)); }
	if (0.5E9 <= 2 * y && y < 1.2E9) { valueResponse = M_HR_pi0 * sigmaHR*(1 - pow((0.5E9) / (2.0*y), 2.0)); }
	if (1.2E9 <= 2 * y) { valueResponse = M_HR_pi0 * sigmaHR*((pow(1.2E9, 2) - pow(0.5E9, 2.0)) / pow(2 * y, 2.0)); }
	return valueResponse;
}


double calculator::calculatorulation()
{
	double funcValue = 0.0;
	double y = 0.0;
	long int countII = 0;
	int count = 0;
	int repetition = 61;           // Anzahl der Proton grid Punkte
	int repetitionPerVec = 3600;   // Q_pi.size()/repetitation

	std::ofstream writingFile;
	writingFile.open("allValues_observe_25.07.18.txt");


	std::ifstream readFile("Ylist_25.07.txt");           // reads file which was created by "calculatorulatesVars()"
	ResultList r = ResultList(60);

	if (!readFile)
	{
		std::cerr << "* Can't open file! *" << std::endl;
	}
	else
	{

		while (readFile >> y >> photonEnergy >> protonEnergy >> pionEnergy)
		{
			if (true)//0.5*145E6 < y)
			{
				funcValue = responseFunc2(y)*photonDistPrep(y);
				valueQ_pi = protonDistPrep(pionEnergy / CHI_pi0LR)*(M_PROTON / pionEnergy)*funcValue;
				Q_pi.push_back(valueQ_pi);
				writingFile << photonEnergy << " " << pionEnergy << " " << protonEnergy << " " << y << " " << responseFunc2(y) << " " << photonDistPrep(y) << std::endl;
			}
			else
			{
				Q_pi.push_back(0.0);
			}
			//std::cout << "reading file and calculatorulate single Q_values." << std::endl;

		}
		std::cout << "finish" << std::endl;


		r.AddValues(Q_pi, repetitionPerVec, repetition);
		std::cout << "filling up vectors. " << std::endl;
	}

	std::cout << "calculatorulating plot values." << std::endl;

	std::ofstream Q_values4Pi;
	Q_values4Pi.open("plotFile4Pions.txt"/*, std::fstream::app*/);

	for (int x = 0; x < r.accumulateQ.size(); x++)
	{
		Q_values4Pi << r.accumulateQ[x] << std::endl;
	}

	std::cout << Q_pi.size() << std::endl;

	writingFile.close();
	Q_values4Pi.close();
	std::cout << "done calculatorulating." << std::endl;
	return 0;
}

