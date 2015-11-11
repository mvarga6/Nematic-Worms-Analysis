/*
This file defines the functor class which launches a
thread to calculate a pair correlation function.

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __SPATIAL_CORRELATOR_HPP__
#define __SPATIAL_CORRELATOR_HPP__

#include <vector>
#include <iostream>

#include "boost\thread.hpp"

#include "litterbox.hpp"
#include "calculator.hpp"
#include "simreplay_types.hpp"
#include "pairs_data.hpp"

namespace nw
{
	template <typename dataType>
	class spatialCorrelator : public calculator
	{
	public:
		spatialCorrelator(){};
		spatialCorrelator(std::vector< xyz<dataType> >& data, pairsData<dataType>& pairData) : target(data), pairs(pairData){};
		spatialCorrelator(std::vector< xyz<dataType> >& data, pairsData<dataType>& pairData, dataType binSize, dataType maxDist) : target(data), pairs(pairData)
		{
			bin_size = binSize; max_dist = maxDist; rho = pairData.getrho();
		};

		//.. calculation specific methods
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

		//.. polymorphic methods
		void start(void)
		{
			std::cout << "G(r) calculation started\n";
			calculationThread = boost::thread(&spatialCorrelator::spatialCorrelatorFunction, this);
		}
		void join(void)
		{
			calculationThread.join();
			std::cout << "G(r) calculation finished\n";
		}

	private:
		//.. members and references
		std::vector<dataType> result;
		std::vector< xyz<dataType> >& target;
		pairsData<dataType>& pairs;

		//.. specifications
		dataType bin_size;
		dataType max_dist;
		dataType rho;

		//.. the meat and potatos
		void spatialCorrelatorFunction();
	};

	//... meat and potatos functions, the calculation at hand
	/*template <typename dataType>
	void spatialCorrelator<dataType>::spatialCorrelatorFunction()
	{
		//.. do the stuff that sets result
		dataType r2_max = max_dist * max_dist;
		unsigned bins = unsigned(ceil(max_dist / bin_size));
		result.resize(bins,0);
		dataType count = 0;
		for (int i = 0; i < target.size(); i++)
		{
			for (int j = i + 1; j < target.size(); j++)
			{
				dataType dx = target[j].x - target[i].x;
				dataType dy = target[j].y - target[i].y;
				dataType dz = target[j].z - target[i].z;
				box.recalcDistance(dx, dy, dz);
				dataType rr = dx*dx + dy*dy + dz*dz;

				if (rr > r2_max || rr < bin_size * bin_size) continue;

				dataType r = sqrt(rr);
				unsigned bin = unsigned(r / bin_size);
				
				if (bin >= result.size() || bin < 0) continue;

				//.. count
				result[bin] += 1;
				
			}
		}

		//.. rescale for area of disk (2*pi*r*dr)
		//	 and total # of contributions (count)
		dataType A = 2.0f * 3.1415f * bin_size * target.size();
		for (int i = 0; i < result.size(); i++)
		{
			result[i] /= (A * (i+1) * bin_size);
		}	
	}*/

	template <typename dataType>
	void spatialCorrelator<dataType>::spatialCorrelatorFunction()
	{
		dataType r2_max = max_dist * max_dist;
		unsigned bins = unsigned(ceil(max_dist / bin_size));
		result.resize(bins, 0);

		//.. loop through all pairs
		dataType count = 0;
		for (int i = 0; i < pairs.size(); i++)
		{
			if (pairs[i].dr >= max_dist) continue;

			count += 1;

			//.. histogram distances
			int bin = int(abs(pairs[i].dr) / bin_size);

			//.. make count
			result.at(bin) += 1;
		}

		//.. rescale for area of disk (2*pi*r*dr)
		//	 and total # of contributions (count)
		//dataType A = /*2.0f * M_PI **/ bin_size * count;
		for (int i = 0; i < result.size(); i++)
		{
			//.. calculate area of ring
			float r1 = (float)i * bin_size;
			float r2 = (float)(i + 1) * bin_size;
			float area = M_PI * (r2*r2 - r1*r1);
			result.at(i) /= area;
			result.at(i) /= rho;
			result.at(i) *= 2.0f;
			result.at(i) /= count;

			std::cout << result.at(i);
			//result.at(i) /= (A * (float(i) + 0.5f) * bin_size);
		}
	}
}

#endif