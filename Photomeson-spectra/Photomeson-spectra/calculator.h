//
// Created by Sergej on 20.06.18.
//
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <numeric>
#include <list>

#ifndef CODE_JUNI_calculator_H
#define CODE_JUNI_calculator_H

#define M_PROTON 938.272E6        //  MeV


class calculator {

public:

	void listProtonEnergy();
	void pionEnergycalculator();
	void sumQ();
	void sumVec();

	double photonDist(double x);
	double photonDistPrep(double x);
	double protonDist(double x);
	double protonDistPrep(double x);
	double responseFunc(double x);
	double responseFunc2(double x);
	double integration();


	double photonDistI(double x);
	double photonDistII(double x);

	double calculatorulation();
	double calculatorulation2();

	void yList();

	void neutrinocalculator();
};


#endif //CODE_JUNI_calculator_H
