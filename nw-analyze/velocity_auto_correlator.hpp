/*
This file defines a velocity auto correlation
calculation functor.

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __VELOCITY_AUTO_CORRELATOR_HPP__
#define __VELOCITY_AUTO_CORRELATOR_HPP__

#include <vector>
#include <random>
#include <iostream>

#include "boost\thread.hpp"

#include "simreplay_types.hpp"
#include "calculator.hpp"
#include "param_set.hpp"

namespace nw
{
	template <typename dataType>
	class velocityAutoCorrelator : public calculator
	{
	public:
		//.. construct with all needed parameters
		velocityAutoCorrelator(parameterSet& paramSet, int binNum, float samplePercent = 1.0f)
		{
			//.. save specifics
			bins = binNum;
			frame_dt = float(paramSet["framerate"].i)*paramSet["dt"].f;
			filled = false;
			result.resize(binNum, 0.0);
			counts_in_bin.resize(binNum, 0);

			//.. calculate random list of targets
			int targs;
			int parts = paramSet["nparticles"].i;
			if ((samplePercent > 0.0f) && (samplePercent < 1.0f))
			{
				targs = int((float)parts * samplePercent);
				for (int i = 0; i < targs; i++)
				{
					//.. choose random target
					int target = rand() % parts;

					//.. check for target in list
					bool found = false;
					for (int j = 0; j < targets.size(); j++)
					{
						if (target == targets[j]) found = true;
					}

					if (found) continue;
					else
					{
						targets.push_back(target);
					}
				}
			}
			//.. list all particles as targets
			else
			{
				targs = parts;
				for (int i = 0; i < targs; i++)
				{
					//.. add all as targets
					targets.push_back(i);
				}
			}

			std::cout << "Done choosing " << targets.size() << " targets\n";
		}

		//.. loads an addition data reference and calculates if ready
		void update(std::vector< xyz<dataType> >& data)
		{
			std::cout << "Reference # " << data_refs.size() << " set with reference\n";

			//.. stays at correct size
			if (filled)
			{
				//.. store data reference
				data_refs.push_back(&data);

				//.. delete 0th element
				data_refs.erase(data_refs.begin());
			}

			//.. first time at correct size
			else 
			{
				//.. store data reference
				data_refs.push_back(&data);

				if (data_refs.size() == bins + 1)
				{
					filled = true;
				}
			}
		}
		void update(std::vector< xyz<dataType> >* data)
		{
			//.. stays at correct size
			std::cout << "Reference # " << data_refs.size() << " set with pointer\n";

			//.. store data reference at end
			data_refs.push_back(data);

			//.. if full
			if (filled)
			{
				std::cout << "data_refs updated\n";

				//.. delete first element
				data_refs.erase(data_refs.begin());
			}

			//.. first time at correct size
			else
			{
				std::cout << "data_refs grows to " << data_refs.size() + 1 << std::endl;

				//.. condition to flag full
				if (data_refs.size() == bins + 1)
				{
					std::cout << "data_refs filled\n";
					filled = true;
				}
			}
		}

		//.. access results and xvalues
		std::vector<dataType> results(void)
		{
			//.. loop through all bins and create average
			for (int i = 0; i < result.size(); i++)
			{
				//.. normallize
				result[i] /= (float)counts_in_bin[i];
				result[i] /= result[0];
			}
			return result;
		}
		std::vector<dataType> xValues(void)
		{
			std::vector<dataType> xvals;
			for (int i = 0; i < bins; i++)
			{
				xvals.push_back(i * frame_dt);
			}
			return xvals;
		}

		//.. polymorphic methods
		void start(void)
		{
			//.. stop if not ready
			if (!filled) return;

			//.. launch thread
			std::cout << "velocity auto-correlation calculation started\n";
			this->calculationThread = boost::thread(&velocityAutoCorrelator::velocityAutoCorrelatorFunction, this);
		}
		void join(void)
		{
			//.. stop if calculate didn't start
			if (!filled) return;

			std::cout << "velocity auto-correlation calculation finished\n";
			this->calculationThread.join();
		}

	private:
		//.. store result
		std::vector<dataType> result;
		std::vector<int> counts_in_bin;

		//.. list of targets' ids
		std::vector<int> targets;

		//.. references to data
		std::vector< std::vector< xyz<dataType> >* > data_refs;
		bool filled;

		//.. specifications
		int bins;
		float frame_dt;
		float sample_percent;

		//.. actual calculation method
		void velocityAutoCorrelatorFunction(void);
	};

	template <typename dataType>
	void velocityAutoCorrelator<dataType>::velocityAutoCorrelatorFunction(void)
	{
		//.. for each target particle
		for (int p = 0; p < targets.size(); p++)
		{
			//.. target id
			int id = targets[p];

			//.. calculate velocities at all time steps for target p
			dataType vx0 = 0.0;
			dataType vy0 = 0.0;
			std::cout << data_refs.size() << std::endl;
			for (int t = 0; t < bins; t++)
			{
				//.. calculate velocity at (t + 1)
				std::cout << "\t" << data_refs[t+1]->size() << std::endl;
				dataType _vx = (data_refs.at(t + 1)->at(id)).x - data_refs.at(t)->at(id).x;
				dataType _vy = (data_refs.at(t + 1)->at(id)).y - data_refs.at(t)->at(id).y;
				_vx /= frame_dt;
				_vy /= frame_dt;

				//.. store if first time
				if (t == 0)
				{
					vx0 = _vx;
					vy0 = _vy;
				}

				//.. add to result the dot products < v(0) * v(dt) >
				result[t] += _vx*vx0 + _vy*vy0;
				counts_in_bin[t]++;
			}
		}
	}
}

#endif