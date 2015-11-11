/*
	This file defines all the available types for
	use with third party program, "Sim Replay".

	M.Varga
	Kent State University
	Chemical Physics Dept.
	2015
*/

#ifndef __SIMREPLAY_TYPES_HPP__
#define __SIMREPLAY_TYPES_HPP__

#include <initializer_list>

namespace nw
{
	class simReplayType
	{
	public:
		virtual void zero(void) = 0;
	};

	template<typename T> class xyzv;
	template<typename T> class xyzc;

	template <typename T>
	class xyz : public simReplayType
	{
	public:
		char t;
		T x, y, z;

		//.. easy init
		xyz(){};
		xyz(const xyz<T>& _xyz) : x(_xyz.x), y(_xyz.y), z(_xyz.z), t(_xyz.t){};
		xyz(const xyzv<T>& _xyzv) : x(_xyzv.x), y(_xyzv.y), z(_xyzv.z), t(_xyzv.t){};
		xyz(std::initializer_list<T> init_list)
		{
			if (init_list.size() == 3)
			{
				x = *(init_list.begin());
				y = *(init_list.begin() + 1);
				z = *(init_list.begin() + 2);
				t = 'A';
			}
		}
		xyz(T _x, T _y, T _z, char _t = 'A') : x(_x), y(_y), z(_z), t(_t){};

		void set(xyz<T>* setter)
		{
			t = setter->t; x = setter->x; y = setter->y; z = setter->z;
		}

		void zero(void)
		{
			t = '\0'; x = y = z = 0;
		}
	};

	template <typename T>
	class xyzc : public simReplayType
	{
	public:
		char t;
		T x, y, z, c;

		//.. easy init
		xyzc(){};
		xyzc(const xyz<T>& _xyz) : x(_xyz.x), y(_xyz.y), z(_xyz.z), t(_xyz.t)
		{
			c = 0;
		}
		xyzc(const xyz<T>& _xyz, T _c) : xyzc(_xyz)
		{
			c = _c;
		}
		xyzc(std::initializer_list<T> init_list)
		{
			if (init_list.size() == 4)
			{
				x = *(init_list.begin());
				y = *(init_list.begin() + 1);
				z = *(init_list.begin() + 2);
				c = *(init_list.begin() + 3);
				t = 'A';
			}
		}
		xyzc(T _x, T _y, T _z, T _c, char _t = 'A') : x(_x), y(_y), z(_z), c(_c), t(_t){};

		void set(xyzc<T>* setter)
		{
			t = setter->t; x = setter->x; y = setter->y; z = setter->z; c = setter->c;
		}

		void zero(void)
		{
			t = '\0'; x = y = z = c = 0;
		}
	};

	template <typename T>
	class xyzv : public simReplayType
	{
	public:
		char t;
		T x, y, z, vx, vy, vz;

		//.. easy init
		xyzv(){};
		xyzv(const xyz<T>& _xyz) : x(_xyz.x), y(_xyz.y), z(_xyz.z), t(_xyz.t)
		{
			vx = vy = vz = 0;
		}
		xyzv(const xyz<T>& _xyz, T _vx, T _vy, T _vz) : xyzv(_xyz)
		{
			vx = _vx; vy = _vy; vz = _vz;
		}
		xyzv(std::initializer_list<T> init_list)
		{
			if (init_list.size() == 6)
			{
				x = *(init_list.begin());
				y = *(init_list.begin() + 1);
				z = *(init_list.begin() + 2);
				vx = *(init_list.begin() + 3);
				vy = *(init_list.begin() + 4);
				vz = *(init_list.begin() + 5);
				t = 'A';
			}
		}
		xyzv(T _x, T _y, T _z, T _vx, T _vy, T _vz, char _t = 'A') : x(_x), y(_y), z(_z), vx(_vx), vy(_vy), vz(_vz), t(_t){};
	
		void set(xyzv<T>* setter)
		{
			t = setter->t; x = setter->x; y = setter->y; z = setter->z; vx = setter->vx; vy = setter->vy; vz = setter->vz;
		}

		void zero(void)
		{
			t = '\0'; x = y = z = vx = vy = vz = 0;
		}
	};
}

#endif