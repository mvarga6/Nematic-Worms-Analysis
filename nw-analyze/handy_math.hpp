/*
This file defines some handy math functions 
useful for all nw namespace stuff

*** Some functions require system
parameters directly from simulation ***

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __HANDY_MATH_HPP__
#define __HANDY_MATH_HPP__

#define _USE_MATH_DEFINES

#include <vector>
#include <math.h>

namespace nw
{
	//.. normalizes the components sent as arguments and
	//	 returns the magnitude as well.
	double normalize(double& u, double& v, double& w)
	{
		double mag = sqrt(u*u + v*v + w*w);

		u /= mag;
		v /= mag;
		w /= mag;

		return mag;
	}
	double normalize(double& u, double& v)
	{
		double mag = sqrt(u*u + v*v);

		u /= mag;
		v /= mag;

		return mag;
	}
	float normalize(float& u, float& v, float& w)
	{
		float mag = sqrtf(u*u + v*v + w*w);

		u /= mag;
		v /= mag;
		w /= mag;

		return mag;
	}
	float normalize(float& u, float& v)
	{
		double mag = sqrt(u*u + v*v);

		u /= mag;
		v /= mag;

		return mag;
	}

	//.. recalculates distances across periodic boundary
	void periodicBC(double& dr, const double& L)
	{
		double Lo2 = L / 2.0;
		if (dr > Lo2) dr -= L;
		if (dr < -Lo2) dr += L;
	}
	void periodicBC(float& dr, const float& L)
	{
		float Lo2 = L / 2.0f;
		if (dr > Lo2) dr -= L;
		if (dr < -Lo2) dr += L;
	}

	//.. gets relative angle
	void angularBC(double& raw_difference)
	{
		if (raw_difference > M_PI) raw_difference -= 2 * M_PI;
		if (raw_difference < -M_PI) raw_difference += 2 * M_PI;
	}
	void angularBC(float& raw_difference)
	{
		if (raw_difference > M_PI) raw_difference -= 2 * M_PI;
		if (raw_difference < -M_PI) raw_difference += 2 * M_PI;
	}

	//.. returns a histogram of vales sent
	template <typename _ty>
	std::vector<_ty> histogramValues(std::vector<_ty>& data, int bins, _ty& width, const _ty& max, const _ty& min)
	{
		//.. empty, resize, and entry amount
		std::vector<_ty> hist;
		hist.resize(bins, 0.0);
		const double entry = 1 / (double)data.size();

		//.. set property
		width = (max - min) / (double)bins;

		//.. add entry to proper bin
		for (int i = 0; i < data.size(); i++)
		{
			int b = int(floor((data.at(i) - min) / width));

			if (b < 0 || b >= bins) continue;
			hist.at(b) += entry;
		}

		return hist;
	}

	template <typename _ty>
	std::vector<_ty> averageData(std::vector< std::vector<_ty> >& data)
	{
		//.. store sizes
		int outerSize = data.size();
		std::vector<int> sizes;
		int maxSize = 0;
		for (int i = 0; i < outerSize; i++)
		{
			int s = data[i].size();
			sizes.push_back(s);
			if (s > maxSize) maxSize = s;
		}

		//.. result is of maxSize
		std::vector<_ty> result;
		result.resize(maxSize, 0.0);

		//.. add and average
		std::vector<double> counts;
		counts.resize(maxSize, 0.0);

		//.. loop over elements
		for (int element = 0; element < maxSize; element++)
		{
			//.. for each vectors
			for (int vec = 0; vec < outerSize; vec++)
			{
				//.. get counts for each element index
				if (element < sizes[vec])
				{
					//.. sum and count
					result[element] += data[vec][element];
					counts[element] += 1.0;
				}
			}
		}

		//.. divide each element by count
		for (int element = 0; element < maxSize; element++)
		{
			result[element] /= counts[element];
		}

		return result;
	}

	template <typename _ty>
	std::vector<_ty> averageData(std::vector< std::vector<_ty>* >& data)
	{
		//.. store sizes
		int outerSize = data.size();
		std::vector<int> sizes;
		int maxSize = 0;
		for (int i = 0; i < outerSize; i++)
		{
			int s = data[i]->size();
			sizes.push_back(s);
			if (s > maxSize) maxSize = s;
		}

		//.. result is of maxSize
		std::vector<_ty> result;
		result.resize(maxSize, 0.0);

		//.. add and average
		std::vector<double> counts;
		counts.resize(maxSize, 0.0);

		//.. loop over elements
		for (int element = 0; element < maxSize; element++)
		{
			//.. for each vectors
			for (int vec = 0; vec < outerSize; vec++)
			{
				//.. get counts for each element index
				if (element < sizes[vec])
				{
					//.. sum and count
					result[element] += data.at(vec)->at(element);
					counts[element] += 1.0;
				}
			}
		}

		//.. divide each element by count
		for (int element = 0; element < maxSize; element++)
		{
			result[element] /= counts[element];
		}

		return result;
	}
}

#endif