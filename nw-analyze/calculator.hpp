/*
This file defines the abstract functor class which 
launches a thread to calculate a desired property
based on the inheriting types.

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __CALCULATOR_HPP__
#define __CALCULATOR_HPP__

#include "boost\thread.hpp"

namespace nw
{
	class calculator
	{
	public:
		virtual void start(void) = 0;
		virtual void join(void) = 0;

	protected:
		boost::thread calculationThread;
	};
}

#endif