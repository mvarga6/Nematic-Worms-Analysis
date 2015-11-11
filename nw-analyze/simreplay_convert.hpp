/*
This file defines the convertions between different
types of Sim Replay types.  

*** Some functions require system 
parameters directly from simulation ***

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __SIMREPLAY_CONVERT_HPP__
#define __SIMREPLAY_CONVERT_HPP__

#include <vector>

#include "simreplay_types.hpp"
//#include "litterbox.hpp"
#include "handy_math.hpp"
#include "param_set.hpp"

using namespace std;

namespace nw
{
	/*
		This funtion returns a vector of xyzv format with
		the v place holder filled with particle orientation 
		dependent on a particles neighbors in the worm.

		*** requires knowledge of particles per worm ***
	*/
	template <typename T>
	vector<xyzv<T>> xyz_to_xyzv(const vector<xyz<T>>& _xyz, parameterSet& params)
	{
		vector<xyzv<T>> result;
		const int partsPerWorm = params["np"].i;
		for (int i = 0; i < _xyz.size(); i++)
		{
			//.. particle number in worm
			unsigned p = i % partsPerWorm;

			//.. for all particles but heads(tails)
			if (p < partsPerWorm - 1)
			{
				//.. calculate distance to ahead neighbor
				T dx = _xyz[i + 1].x - _xyz[i].x;
				T dy = _xyz[i + 1].y - _xyz[i].y;

				//.. BC's and unit vector
				// BOUNDARY CONDITION FUNCTION GOES HERE !!!!
				nw::normalize(dx, dy);

				//.. add to result vector
				xyzv<T> tmp(_xyz[i], dx, dy, 0.0);
				result.push_back(tmp);
			}
			//.. the heads(tails) take on the previous angle
			else
			{
				xyzv<T> tmp(_xyz[i], result[i - 1].vx, result[i - 1].vy, result[i - 1].vz);
				result.push_back(tmp);
			}
		}

		//.. return created vector
		return result;
	}

	template <typename T>
	vector<xyzc<T>> xyz_to_xyzc(const vector<xyz<T>>& _xyz, const vector<T>& cValues)
	{
		vector<xyzc<T>> result;

		//.. execute if sizes match
		if (_xyz.size() == cValues.size())
		{
			for (int i = 0; i < _xyz.size(); i++)
			{
				xyzc<T> tmp(_xyz[i], cValues[i]);
				result.push_back(tmp);
			}
		}
		//.. otherwise send NULLed array back
		else
		{
			result.resize(_xyz.size(), NULL);
		}
		
		return result;
	}

	template <typename T>
	void randomShrink(vector<xyz<T>>& _xyz, int numberToRemove)
	{
		for (int i = 0; i < numberToRemove)
		{
			int rid = rand() % _xyz.size();
			_xyz.erase(rid);
		}
	}
}
#endif