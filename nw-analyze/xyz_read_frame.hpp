/*
	This file defines all the available ways to read
	a SimReplay *.xyz file.

	M.Varga
	Kent State University
	Chemical Physics Dept.
	2015
*/

#ifndef __XYZ_READ_HPP__
#define __XYZ_READ_HPP__

#include <fstream>
#include <vector>
#include "simreplay_types.hpp"

using namespace std;

namespace nw
{
	/*
		Returns false if function reaches end of file flag or
		send particle number doesn't agree with file number.
		All necessary parameters sent as arguements.
		Loads the arrays of known size with particle positions
		excluding the corner particles for one frame.

		Safe Method! i.e. ifstream safetly checks build in!

		--------- USAGE ----------
		T:				type of data (wanted) in containers
		fin:			pre-open in file stream
		c,x,y,z:		containers to be loaded with positions and types
		partsInFile:	number of particles defined in file (excluding corners)
		preAllocated:	determines if arrays need memory allocation
		cornersNumber:	to know how many particle positions to trash at end of frame (default = 4)
	*/
	template <typename dataType>
	bool xyz_read_frame(ifstream& fileStream, char* c, dataType* x, dataType* y, dataType* z, const unsigned& particleNumber, bool preAllocated = true, unsigned cornersNumber = 4)
	{
		//.. safety and null 
		if (!fileStream.is_open())
		{
			x = y = z = nullptr;
			c = nullptr;
			return false;
		}

		//.. tmp variables
		unsigned parts, np;
		float hx, hy, k2;
		dataType Ttrash;
		char ctrash;

		//.. check for eof
		if (fileStream.eof()) return false;

		//.. load particle number
		fileStream >> parts;

		//.. check for eof
		if (fileStream.eof()) return false;

		//.. read in system parameters (needs a safe function)
		fileStream >> np >> k2 >> hx >> hy;

		//.. stop if sizes dont agree
		if (parts != particleNumber + 4) return false;

		//.. allocated new memory if needed
		if (!preAllocated)
		{
			c = new char[particleNumber];
			x = new dataType[particleNumber];
			y = new dataType[particleNumber];
			z = new dataType[particleNumber];
		}

		//..  load all data
		for (int i = 0; i < particleNumber; i++)
		{
			if (!fileStream.eof())
			{
				fileStream >> c[i] >> x[i] >> y[i] >> z[i];
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

		//.. successful read
		return true;
	}

	/*
		Returns false if function reaches end of file flag.
		Parameters inferred from file.
		Allocates Memory based on frame in file and returns 
		the size of allocated ptr arrays in particleNumber param. 
		Loads the arrays with particle positions excluding the 
		corner particles for single frame of file.
	
		Unsafe Method! i.e. NO ifstream safetly checks build in!

		--------- USAGE ----------
		T:				type of data (wanted) in containers
		fin:			pre-open in file stream
		c,x,y,z:		containers to be loaded with positions and types
		particleNumber: sets the particle number by reference
		cornersNumber:	to know how many particle positions to trash at end of frame (default = 4)
	*/
	template <typename dataType>
	bool xyz_read_frame(ifstream& fileStream, char* c, dataType* x, dataType* y, dataType* z, unsigned& particleNumber, unsigned cornersNumber = 4)
	{
		//.. tmp variables
		unsigned parts, np, pcount = 0;
		float hx, hy, k2;
		dataType Ttrash;
		char ctrash;

		//.. check end of file
		if (fileStream.eof()) return false;

		//.. read in particle number
		fileStream >> parts;

		//.. check end of file
		if (fileStream.eof()) return false;

		//.. read in system parameters  (needs a safe function)
		fileStream >> np >> k2 >> hx >> hy;

		//.. correct for corners
		parts -= cornersNumber;

		//.. allocated new memory always
		c = new char[parts];
		x = new dataType[parts];
		y = new dataType[parts];
		z = new dataType[parts];

		//..  load all data
		for (int i = 0; i < parts; i++)
		{
			if (!fileStream.eof())
			{
				fileStream >> c[i] >> x[i] >> y[i] >> z[i];
				pcount++;
			}
			else
			{
				particleNumber = pcount;
				return false;
			}
		}

		//.. trash the corner particles
		for (int i = 0; i < cornersNumber; i++)
		{
			if (!fileStream.eof())
			{
				fileStream >> ctrash >> Ttrash >> Ttrash >> Ttrash;
			}
			else
			{
				particleNumber = pcount;
				return false;
			}
		}

		//.. assign particle number based on frame
		particleNumber = parts;

		//.. finished without eof
		return true;
	}

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
		c,x,y,z:		containers to be loaded with positions
		cornersNumber:	to know how many particle positions to trash at end of frame (default = 4)
	*/
	template <typename dataType>
	bool xyz_read_frame(ifstream& fileStream, vector<char>& c, vector<dataType>& x, vector<dataType>& y, vector<dataType>& z, unsigned cornersNumber = 4)
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
		c.clear();
		x.clear();
		y.clear();
		z.clear();

		//..  load all data
		for (int i = 0; i < parts; i++)
		{
			if (!fileStream.eof())
			{
				dataType _x, _y, _z;
				char _c;
				fileStream >> _c >> _x >> _y >> _z;
				c.push_back(_c); x.push_back(_x); y.push_back(_y); z.push_back(_z);
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
	bool xyz_read_frame(ifstream& fileStream, vector<xyz<dataType>>& data, unsigned cornersNumber = 4)
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
				dataType _x, _y, _z;
				char _t;
				fileStream >> _t >> _x >> _y >> _z;
				xyz<dataType> _r(_x, _y, _z, _t);
				data.push_back(_r);
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