/*
This file defines a container for a linked list 
designed to be used for faster neighbor finding.

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __LINKED_LIST_HPP__
#define __LINKED_LIST_HPP__

#include <vector>
#include <math.h>
#include <iostream>

#include "simreplay_types.hpp"
#include "param_set.hpp"
#include "pairs_data.hpp"

namespace nw
{
	template <typename ty>
	class linkedList
	{
		template <typename dataType>
		friend class pairsData;

	public:

		//.. create with data and system parameters
		linkedList(std::vector< xyz<ty> >& data, parameterSet& params, float cellSize)
		{
			xbox = params["xbox"].f;
			ybox = params["ybox"].f;
			np = params["np"].i;
			dcell = cellSize;
			rho = (params["nworms"].i * ((np - 1)*params["l1"].f + params["sigma"].f * (pow(2.0f, 1.0f / 6.0f)))) / (xbox * ybox);

			std::cout << rho << std::endl;

			nxcell = int(ceil(xbox / dcell));
			nycell = int(ceil(ybox / dcell));
			ncells = nxcell * nycell;

			this->generate(data);
		}
		
		//.. set properties
		void setProperties(parameterSet& params, float cellSize)
		{
			xbox = params["xbox"].f;
			ybox = params["ybox"].f;
			np = params["np"].i;
			dcell = cellSize;

			nxcell = int(ceil(xbox / dcell));
			nycell = int(ceil(ybox / dcell));
			ncells = nxcell * nycell;
		}

		//.. resets list with same parameters
		void reset(std::vector< xyz<ty> >& data)
		{
			//.. store new data reference
			this->generate(data);
		}
		
		friend class pairsData<ty>;

		void printHeads(void)
		{
			for (std::vector<int>::iterator it = heads.begin(); it != heads.end(); ++it)
				std::cout << it - heads.begin() << ": " << *it << std::endl;
		}
		void printPointers(void){
			for (std::vector<int>::iterator it = pointers.begin(); it != pointers.end(); ++it)
				std::cout << it - pointers.begin() << ": " << *it << std::endl;
		}
		int  countHeads(void)
		{
			int count = 0;
			for (std::vector<int>::iterator it = pointers.begin(); it != pointers.end(); ++it)
			{
				if (*it == -1) count++;
			}
			return count;
		}

	private:
		//.. lists
		std::vector<int> heads;    // sized by specs
		std::vector<int> pointers; // sized by data
		//std::vector< xyz<ty> >& data_ref;

		//.. specifications
		int nxcell, nycell, ncells;

		//.. required system parameters
		ty xbox, ybox, dcell, rho;
		int np;

		//.. generates lists
		void generate(std::vector< xyz<ty> >& data_ref)
		{
			std::cout << "Generating linked list\t";

			//.. allocate
			heads.resize(ncells, -1);
			pointers.resize(data_ref.size(), -1);

			//.. loop through all particles
			for (int i = 0; i < data_ref.size(); i++)
			{
				//.. find targets cell
				int icell = int(floor(data_ref[i].x / dcell));
				int jcell = int(floor(data_ref[i].y / dcell));
				int scell = jcell*nxcell + icell;

				//.. set point to previous head of that cell
				pointers[i] = heads[scell];

				//.. set new head of that cell to that particle
				heads[scell] = i;
			}

			std::cout << "finished\n";
		}
	};
}

#endif