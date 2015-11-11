/*
This file defines a calculation for making a 
histogram of particle orientations.

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __ORIENTATION_HIST_HPP__
#define __ORIENTATION_HIST_HPP__

#include <vector>
#include <math.h>

#include "boost\thread.hpp"

#include "handy_math.hpp"
#include "calculator.hpp"

namespace nw
{
	template <typename dataType>
	class orientationHistogram : public calculator
	{
	public:

		//.. create with all needed things
		orientationHistogram(std::vector< xyzv<dataType> >& data, int numberOfBins) : target(data), bins(numberOfBins){};
		
		//.. getting access to calculation
		std::vector<dataType> results(void)
		{
			return hist;
		}
		std::vector<dataType> xValues(void)
		{
			std::vector<dataType> xvals;
			for (int i = 0; i < hist.size(); i++)
			{
				xvals.push_back(i*bin_width + min);
			}
			return xvals;
		}

		//.. polymorphic methods
		void start(void)
		{
			std::cout << "Histogram being produced\n";
			calculationThread = boost::thread(&orientationHistogram::orientationHistogramFunction, this);
		}
		void join(void)
		{
			calculationThread.join();
			std::cout << "Histogram finished\n";
		}

	private:
		//.. members and references
		std::vector< xyzv<dataType> >& target;
		std::vector<dataType> hist;

		//.. specifications
		int bins;
		dataType max, min;
		dataType bin_width;

		//.. calculation function
		void orientationHistogramFunction();
	};

	template <typename dataType>
	void orientationHistogram<dataType>::orientationHistogramFunction()
	{
		//.. produce thetas and find max and min
		max = 0.0; min = 0.0;
		std::vector<dataType> theta;
		for (int i = 0; i < target.size(); i++)
		{	
			dataType t = atan2(target[i].vy, target[i].vx);
			theta.push_back(t);
			if (t < min) min = t;
			if (t > max) max = t;
		}

		//.. make histogram (sets bin_width)
		hist = nw::histogramValues<dataType>(theta, bins, bin_width, max, min);
	}
}

#endif