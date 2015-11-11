/*
This file defines the functor class which launches a
thread to calculate an orientation correlation.

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __ORIENTATION_CORRELATOR_HPP__
#define __ORIENTATION_CORRELATOR_HPP__

#include <vector>
#include <iostream>

#include "boost\thread.hpp"
#include "simreplay_types.hpp"
#include "pairs_data.hpp"
#include "calculator.hpp"

namespace nw
{
	template <typename dataType>
	class orientationCorrelator : public calculator
	{
	public:
		orientationCorrelator(){};
		orientationCorrelator(std::vector< xyzv<dataType> >& data, pairsData<dataType>& pairData) : target(data), pairs(pairData){};
		orientationCorrelator(std::vector< xyzv<dataType> >& data, pairsData<dataType>& pairData, dataType binSize, dataType maxDist) : target(data), pairs(pairData)
		{
			this->setProperties(binSize, maxDist);
		};

		void setProperties(dataType binSize, dataType maxDist)
		{
			bin_size = binSize; max_dist = maxDist;
		}
		std::vector<dataType> results(void)
		{
			return result;
		}
		std::vector<dataType> xValues(void)
		{
			std::vector<dataType> xvals;
			for (int i = 0; i < result.size(); i++)
			{
				xvals.push_back(i*bin_size);
			}
			return xvals;
		}

		void start(void)
		{
			std::cout << "OCF(r) calculation started\n";
			calculationThread = boost::thread(&orientationCorrelator::orientationCorrelatorFunction, this);
		}
		void join(void)
		{
			calculationThread.join();
			std::cout << "OCF(r) calculation finshed\n";
		}

	private:
		//.. members and references
		std::vector<dataType> result;
		std::vector< xyzv<dataType> >& target;
		pairsData<dataType>& pairs;

		//.. specifications
		dataType bin_size;
		dataType max_dist;

		//.. meat and potatos
		void orientationCorrelatorFunction();
	};

	template <typename dataType>
	void orientationCorrelator<dataType>::orientationCorrelatorFunction()
	{
		dataType r2_max = max_dist * max_dist;
		unsigned bins = unsigned(ceil(max_dist / bin_size));
		result.resize(bins, 0);

		//.. loop through pairs
		dataType count = 0;
		for (int i = 0; i < pairs.size(); i++)
		{
			dataType r = pairs[i].dr;
			if (r >= max_dist) continue;

			count += 1;

			//.. find bin
			int bin = int(abs(r) / bin_size);

			//.. sum dot products
			result.at(bin) += abs(pairs[i].dot());
		}

		//.. rescale for area of binned disks
		dataType A = /*2.0f * M_PI **/ bin_size * count;
		for (int i = 0; i < result.size(); i++)
		{
			result.at(i) /= (A * ((float)i + 0.5f) * bin_size);
		}
	}
}

#endif