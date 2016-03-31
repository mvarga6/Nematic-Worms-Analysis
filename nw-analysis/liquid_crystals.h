// Nematic Worms Analysis
// 2.4.16
// liquid crystal generics
// Mike Varga
#ifndef __LIQUID_CRYSTALS_H__
#define __LIQUID_CRYSTALS_H__
// -------------------------------
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
// -------------------------------
class n {
public:

	inline n(){
		nx_ = 0;
		ny_ = 0;
		nz_ = 0;
		mag_ = 0;
	}
	inline n(float x, float y, float z = 0){
		nx_ = x;
		ny_ = y;
		nz_ = z;
		mag_ = std::sqrt(x*x + y*y + z*z);
	}
	inline void set(float x, float y, float z = 0){
		nx_ = x;
		ny_ = y;
		nz_ = z;
		mag_ = std::sqrt(x*x + y*y + z*z);
	}
	inline float mag(){ return mag_; }
	inline void normalize(){
		mag_ = sqrt(nx_*nx_ + ny_*ny_ + nz_*nz_);
		if (mag_ == 0.0f) return;
		nx_ /= mag_;
		ny_ /= mag_;
		nz_ /= mag_;
	}
	inline float& getnx(){ return nx_; }
	inline float& getny(){ return ny_; }
	inline float& nx(){ return nx_; }
	inline float& ny(){ return ny_; }
	inline float& nz(){ return nz_; }
	float& operator[](int d) {
		if (d == 0) return nx_;
		else if (d == 1) return ny_;
		else if (d == 2) return nz_;
		else return mag_;
	}
	inline float operator*(const n& rhs){
		return (this->nx_*rhs.nx_ + this->ny_*rhs.ny_ + this->nz_*rhs.nz_);
	}
	inline float cos_angle_between(n& n2){
		return ((*this) * n2) / (this->mag() * n2.mag());
	}
	inline float angle_between(n& n2){
		float ang = acos(this->cos_angle_between(n2));
		if (ang > M_PI) ang -= 2.0f * M_PI;
		if (ang < -M_PI) ang += 2.0f * M_PI;
		return ang;
	}
	
private:
	float nx_;
	float ny_;
	float nz_;
	float mag_;
};

typedef std::vector<n> ngroup;
typedef ngroup::iterator ngroup_iter;
float calculate_S(ngroup& u, bool true3D_S = false) {
	if (u.empty()) return 0; // handle empty list
	for (auto &it : u) it.normalize(); // normalize all
	n dir(0.0f, 0.0f, 0.0f);
	for (auto &it : u) { // formulate director
		dir.nx() += it.nx();
		dir.ny() += it.ny();
		dir.nz() += it.nz();
	}
	dir.normalize(); // make unit vector
	if (dir.mag() <= 0.0f) return 0;
	float cos2ave = 0.0f, cossinave = 0.0f; // accumulated for averaging
	float costh, sinth, cos2th; // middlemen variables
	for (auto &it : u) { // calculate <cos^2(th)>
		costh = dir.cos_angle_between(it);
		cos2th = costh*costh;
		sinth = sqrt(1 - cos2th);
		cos2ave += cos2th;
		cossinave += costh*sinth;
	}
	cos2ave /= float(u.size()); // make numeric averages
	cossinave /= float(u.size());
	if (true3D_S)
		return (0.5f * (3 * cos2ave - 1)); // return true 3D order
	return (2.0f * (sqrt((cos2ave - 0.5f)*(cos2ave - 0.5f) + cossinave*cossinave)));
}

#endif