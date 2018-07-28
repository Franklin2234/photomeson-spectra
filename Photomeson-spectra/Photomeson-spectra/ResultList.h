#include <list>
#include <vector>

#pragma once
class ResultList
{
private:
	int size;
	
public:
	std::vector<std::vector<double>> list;
	std::vector<double> accumulateQ;

	ResultList(int size);
	~ResultList();
	
	void AddValues(std::vector<double> vec, int repetitionPerVec, int repetition);
	void Accumulate();
};