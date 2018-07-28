#include "ResultList.h"
#include "calculator.h"

ResultList::ResultList(int size)
{
	this->size = size;
	for (int x = 0; x < size; x++)
	{
		list.push_back(std::vector<double>());
	}
}


ResultList::~ResultList()
{
}

void ResultList::AddValues(std::vector<double> vec, int repetitionPerVec, int repetition)
{
	for (int x = 0; x < repetitionPerVec; x++)
	{
		for (int y = 0; y < size; y++)
		{
			list[y].push_back(vec[y + x * repetition]);
		}
	}
	Accumulate();
}

void ResultList::Accumulate()
{
	for (int y = 0; y < size; y++)
	{
		accumulateQ.push_back(std::accumulate(list[y].begin(), list[y].end(), 0.0));
	}
}
