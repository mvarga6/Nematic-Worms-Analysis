/*
This file defines a storage class with all info
about pair-wise data.

i.e. Neighbors lists which include relative orientation,
	 distance (if within defined cutoff), and relative
	 velocities.

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __PAIRS_DATA_HPP__
#define __PAIRS_DATA_HPP__

#include <math.h>
#include <vector>
#include <random>
#include <stdlib.h>

#include "simreplay_types.hpp"
//#include "litterbox.hpp"
#include "linked_list.hpp"
#include "param_set.hpp"
#include "handy_math.hpp"

namespace nw
{
	//static const int ddx[9] = {0,-1,0,1,1,1,0,-1,-1};
	//static const int ddy[9] = {0,1,1,1,0,-1,-1,-1,0};

	//.. produces vectors of indices for neighboring boxes
	//	 return size for looping
	int setupNeighborsLoop(std::vector<int>& ddx, std::vector<int>& ddy, double dcell, double Rmax)
	{
		//.. clean arrays
		ddx.empty();
		ddy.empty();

		//.. cellular radius
		int M = int(ceil(Rmax / dcell));
		int MM = M*M;

		//.. loop over i index
		for (int i = -M; i <= M; i++)
		{
			//.. calculate number of j's at that i
			int jmax = int(ceil(sqrt(MM - i*i)));

			//.. loop j+1 times for adding
			for (int j = 0; j <= jmax; j++)
			{
				ddx.push_back(i);
				ddy.push_back(j);
			}
		}

		//.. return sizes
		return ddx.size();
	}

	//.. info for one pair
	template <typename _ty>
	struct pair
	{
		_ty dr;  //.. distance between
		_ty phi; //.. relative orientation
		int a;   //.. index one
		int b;   //.. index two

		pair(_ty _dr, _ty _phi, int _a, int _b) : dr(_dr), phi(_phi), a(_a), b(_b){};
		_ty dot(void)
		{
			return cos(phi);
		}
	};

	//.. container for all pair-wise information
	template <typename dataType>
	class pairsData
	{
	public:
		//.. initialize with everything then run (slower searching)
		pairsData(std::vector< xyzv<dataType> >& data, double searchCutoff, parameterSet& paramSet, bool intraWorm = false)
		{
			//.. specifications
			cutoff = searchCutoff;
			intra = intraWorm;

			//.. saving needed parameters
			xbox = paramSet["xbox"].f;
			ybox = paramSet["ybox"].f;
			np = paramSet["np"].i;

			//.. run with all info accumulated
			this->run(data, searchCutoff*searchCutoff);
		}
		
		//.. initialize with Linked List (faster searching)
		pairsData(std::vector< xyzv<dataType> >& data, double searchCutoff, linkedList<dataType>& list, int skipEvery = -1, bool intraWorm = false)
		{
			cutoff = searchCutoff;
			intra = intraWorm;
			xbox = list.xbox;
			ybox = list.ybox;
			np = list.np;
			rho = list.rho;

			this->runUsingList(data, cutoff*cutoff, list, skipEvery);
		}

		//.. change the parameterSet and criteria
		void setProperties(double searchCutoff, parameterSet& paramSet, bool intraWorm = false)
		{
			//.. specifications
			cutoff = searchCutoff;
			intra = intraWorm;

			//.. saving needed parameters
			xbox = paramSet["xbox"].f;
			ybox = paramSet["ybox"].f;
			np = paramSet["np"].i;
		}
		
		//.. rerun the pairs finding function
		void reset(std::vector< xyzv<dataType> >& data)
		{
			//.. eliminate current container
			pairs.empty();
			this->run(data, cutoff*cutoff);
		}

		//.. rerun the pairs finding function using a linked list
		void reset(std::vector< xyzv<dataType> >& data, linkedList<dataType>& list)
		{
			//.. eliminate current container
			pairs.empty();
			this->runUsingList(data, cutoff*cutoff, list);
		}
		
		dataType getrho(){ return rho; }

		//.. returns number of pairs
		size_t size(void)
		{
			return pairs.size();
		}

		//.. overloading for easy access
		pair<dataType>& operator[](int i)
		{
			return pairs.at(i);
		}

	private:
		//.. pairs data
		std::vector< pair<dataType> > pairs;

		//..  specifications
		double cutoff;
		bool intra;

		//.. encapsolation of needed parameters
		dataType xbox, ybox, rho;
		int np;

		//.. generates the list by searching
		void run(std::vector< xyzv<dataType> >& data_ref, double r2cut);

		//.. generates the list using linkedList method
		void runUsingList(std::vector< xyzv<dataType> >& data_ref, double r2cut, linkedList<dataType>& list, int skipEvery = -1);
	};

	template <typename dataType>
	void pairsData<dataType>::run(std::vector< xyzv<dataType> >& data_ref, double r2cut)
	{
		//.. loop over all pairs
		std::cout << "Searching for all pairs.\n";
		for (int i = 0; i < data_ref.size() - 1; i++)
		{
			for (int j = i + 1; j < data_ref.size(); j++)
			{
				//.. reject intra worm pairs if desired
				if (!intra)
				{
					int wi = i / np;
					int wj = j / np;
					if (wi == wj) continue;
				}

				//.. distances
				dataType dx = data_ref[j].x - data_ref[i].x;
				dataType dy = data_ref[j].y - data_ref[i].y;
				
				//.. boundary conditions
				nw::periodicBC(dx, xbox);
				nw::periodicBC(dy, ybox);

				//.. dist square and reject condition
				dataType rr = dx*dx + dy*dy;
				if (rr > r2cut) continue;

				//.. relative orientation (using smallest in -PI -> PI)
				dataType dt = atan2(data_ref[j].vy, data_ref[j].vx) - atan2(data_ref[i].vy, data_ref[i].vx);
				angularBC(dt);

				//.. add to list
				pair<dataType> tmp(sqrt(rr), dt, i, j);
				pairs.push_back(tmp);
			}	
		}
		std::cout << "There are " << pairs.size() << " pairs.\n";
	}

	template <typename dataType>
	void pairsData<dataType>::runUsingList(std::vector< xyzv<dataType> >& data_ref, double r2cut, linkedList<dataType>& list, int skipEvery = -1)
	{
		std::cout << "Generating pairs data\n";

		int skip;
		if (skipEvery <= 1)
			skip = data_ref.size();
		else
			skip = skipEvery;


		//.. neighboring cells list
		std::vector<int> ddx, ddy;
		int dirmax = setupNeighborsLoop(ddx, ddy, list.dcell, cutoff);

		//.. print list for user
		for (int i = 0; i < dirmax; i++)
		{
			std::cout << ddx[i] << ", " << ddy[i] << std::endl;
		}

		const int imax = list.nxcell;
		const int jmax = list.nycell;
		for (int ic = 0; ic < imax; ic++)
		{
			for (int jc = 0; jc < jmax; jc++)
			{
				//.. scalar address and stop if cell is empty
				int scell = jc*imax + ic;
				std::cout << scell << " ";
				if (list.heads[scell] == -1) continue;

				//.. loop about neighboring cells
				//for (int dir = 0; dir < 5; dir++)
				for (int dir = 0; dir < dirmax; dir++)
				{
					//.. periodic boundary conditions
					int icnab = (ic + ddx[dir]) % imax;
					int jcnab = (jc + ddy[dir]) % jmax;
					if (icnab < 0) icnab += imax;
					if (jcnab < 0) jcnab += jmax;

					//.. scalar address of neighboring cell and stop if empty
					int scnab = jcnab*imax + icnab;
					if (list.heads[scnab] == -1) continue;

					//.. get id of head of scell
					int ii = list.heads[scell];

					//.. loop through all target cell parts
					while (ii > -1)
					{
						if (rand() % skip == 0)
						{
							//.. setup for next jj iteraction
							ii = list.pointers[ii];
							continue;
						}

						//.. get head of neighboring cell
						int jj = list.heads[scnab];

						//.. loop through all neighboring cell parts
						while (jj > -1)
						{
							//.. when excluding intraworms
							if (!intra)
							{
								int iworm = ii / np;
								int jworm = jj / np;
								if (iworm == jworm)
								{
									//.. setup for next jj iteraction
									jj = list.pointers[jj];
									continue;
								}
							}

							//.. distances
							dataType dx = data_ref[jj].x - data_ref[ii].x;
							dataType dy = data_ref[jj].y - data_ref[ii].y;

							//.. boundary conditions
							periodicBC(dx, xbox);
							periodicBC(dy, ybox);

							//.. dist square and reject condition
							dataType rr = dx*dx + dy*dy;
							if (rr <= r2cut)
							{
								//.. relative orientation (using smallest in -PI -> PI)
								dataType dt = atan2(data_ref[jj].vy, data_ref[jj].vx) - atan2(data_ref[ii].vy, data_ref[ii].vx);
								angularBC(dt);

								//.. add to list
								pair<dataType> tmp(sqrt(rr), dt, ii, jj);
								pairs.push_back(tmp);
							}
							//.. move to next particle in scnab
							jj = list.pointers[jj];
						}
						//.. move to next particle in scell
						ii = list.pointers[ii];
					}

				}
			}
		}

		std::cout << "finished\n";
	}
}

#endif