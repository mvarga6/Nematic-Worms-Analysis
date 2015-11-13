/*
This file defines the container for all the simulation
parameter configurations.  Loaded from files (*_cfg.dat)

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __PARAM_SET_HPP__
#define __PARAM_SET_HPP__

#include <string>
#include <fstream>
#include <map>

namespace nw
{

	union parameter
	{
		int i;
		float f;
		bool b;
	};


	class parameterSet
	{
	public:

		//.. all parameters storage
		std::map<std::string, parameter> set;

		//.. Simulation key name
		std::string simkey;

		parameterSet();
		parameterSet(std::string);

		//.. load parameters
		bool LoadWithStream(std::ifstream*);
		bool LoadWithName(std::string);

		//.. print to console
		void Show(void);
		parameter operator[](std::string key) const;

		~parameterSet();
	};
}
#endif
