/*
This file defines the litterbox class.  It holds all
the information about the simulation box.

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __LITTERBOX_HPP__
#define __LITTERBOX_HPP__

#include <initializer_list>
#include "simreplay_types.hpp"

namespace nw
{
	/* 
		The abstract class defining a type of system box
	*/
	template <typename dataType>
	class litterBox
	{
	protected:
		dataType Lx, Ly, Lz;
	
	public:	
		virtual void executeBC(xyz<dataType>& _xyz) = 0;
		virtual void recalcDistance(dataType& dx, dataType& dy, dataType& dz) = 0;

		litterBox() : Lx(0), Ly(0), Lz(0){};
		litterBox(dataType _Lx, dataType _Ly, dataType _Lz) : Lx(_Lx), Ly(_Ly), Lz(_Lz){};
	};

	/*
		A box with periodic boundary conditions
	*/
	template <typename dataType>
	class periodicBox : public litterBox <dataType>
	{
	private:
		dataType Lxo2, Lyo2, Lzo2;

	public:
		periodicBox() : litterBox(){};
		periodicBox(dataType _Lx, dataType _Ly, dataType _Lz) : litterBox(_Lx, _Ly, _Lz), Lxo2(_Lx / 2.0), Lyo2(_Ly / 2.0), Lzo2(_Lz / 2.0){};
	
		void executeBC(xyz<dataType>& _xyz)
		{
			if (_xyz.x > Lx) _xyz.x -= Lx;
			if (_xyz.y > Ly) _xyz.y -= Ly;
			if (_xyz.z > Lz) _xyz.z -= Lz;

			if (_xyz.x < 0.0) _xyz.x += Lx;
			if (_xyz.y < 0.0) _xyz.y += Ly;
			if (_xyz.z < 0.0) _xyz.z += Lz;
		}
		void recalcDistance(dataType& dx, dataType& dy, dataType& dz)
		{
			if (dx > Lxo2) dx -= Lx;
			if (dy > Lyo2) dy -= Ly;
			if (dz > Lzo2) dz -= Lz;

			if (dx < -Lxo2) dx += Lx;
			if (dy < -Lyo2) dy += Ly;
			if (dz < -Lzo2) dz += Lz;
		}
		
	};

}

#endif