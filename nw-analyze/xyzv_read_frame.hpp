/*
This file defines all the available ways to read
a SimReplay *.xyz file.

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __XYZV_READ_FRAME_HPP__
#define __XYZV_READ_FRAME_HPP__

#include <fstream>
#include <vector>
#include "simreplay_types.hpp"

namespace nw
{
	/*
	Returns false if eof flag reached or file isn't open.
	Parameters inferred from frame in file.
	Erases current vector containers and loads them with
	the particle positions and types in file if doesn't false
	first. Resizes the vectors based on parameters defined in
	file.  Doesn't return size because containers contain it.

	Safe Method! i.e. ifstream safetly checks build in!

	--------- USAGE ----------
	T:				type of data (wanted) in containers
	fin:			pre-open in file stream
	data:		containers to be loaded with positions
	cornersNumber:	to know how many particle positions to trash at end of frame (default = 4)
	*/
	template <typename dataType>
	bool xyzv_read_frame(ifstream& fileStream, vector< xyzv<dataType> >& data, unsigned cornersNumber = 4)
	{
		if (!fileStream.is_open())
		{
			return false;
		}

		//.. tmp variables
		unsigned parts, np;
		float hx, hy, k2;
		dataType Ttrash;
		char ctrash;

		//.. check for eof
		if (fileStream.eof()) return false;

		//.. load particles number
		fileStream >> parts;

		//.. check for eof
		if (fileStream.eof()) return false;

		//.. read in comment line (needs a safe function)
		fileStream >> np >> k2 >> hx >> hy;

		//.. correct for corners
		parts -= cornersNumber;

		//.. erase current stuff in containers
		data.clear();

		//..  load all data
		for (int i = 0; i < parts; i++)
		{
			if (!fileStream.eof())
			{
				dataType _x, _y, _z, _vx, _vy, _vz;
				char _t;
				fileStream >> _t >> _x >> _y >> _z >> _vx >> _vy >> _vz;
				xyzv<dataType> _rv(_x, _y, _z, _vx, _vy, _vz, _t);
				data.push_back(_rv);
			}
			else return false;
		}

		//.. trash the corner particles
		for (int i = 0; i < cornersNumber; i++)
		{
			if (!fileStream.eof())
			{
				fileStream >> ctrash >> Ttrash >> Ttrash >> Ttrash;
			}
			else return false;
		}

		//.. successful reaad
		return true;
	}
}

#endif