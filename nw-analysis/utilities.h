#pragma once
// Nematic Worms Analysis
// 2.9.16
// Utility Functions
// Mike Varga
// ------------------------------
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <algorithm>
// ------------------------------
namespace util {
	
	//.. utilities for reading and writing *.xyz files
	namespace simReplay {

		//.. definitions
		typedef std::string str;
		typedef std::stringstream sstrm;
		typedef const str cstr;
		class properties;

		//.. propotype of reading lines
		bool readLine(std::ifstream&, char&, float&, float&, float&);
		bool readLine(std::ifstream&, char&, float&, float&, float&, float&);
		bool readLine(std::ifstream&, char&, float&, float&, float&, float&, float&, float&);
		bool getProperties(str&, properties&);

		//.. types of simreplay files
		enum class type {
			xyz,
			xyzc,
			xyzv,
			no_type
		};

		//.. properties of simreplay files
		class properties {
			bool set_; // state of properties setting
			int atoms_; // number of atoms including corners
			int aggsize_; // inferred aggregate size (from atom types)
			int cornrs_; // number of corners (default 4)
			str comment_; // comment line
			type type_; // file type xyz, xyzc, xyzv
			float xmin_, xmax_; // min-max for dims
			float ymin_, ymax_;
			float zmin_, zmax_;
		public:
			properties(int atms = 0, str cmt = "", type typ = type::no_type)
				: atoms_(atms), comment_(cmt), type_(typ), cornrs_(4), set_(false),
				xmin_(0), xmax_(0), ymin_(0), ymax_(0), zmin_(0), zmax_(0), aggsize_(1){};
			~properties(){};

			//.. set properties from file
			bool getFrom(str& file) { return getProperties(file, *this); }
			
			bool is_set() { return this->set_; }
			void set() { this->set_ = true; }
			int& atoms(){ return this->atoms_; }
			void atoms(int a) { this->atoms_ = a; };
			int& aggsize(){ return this->aggsize_; }
			void aggsize(int a) { this->aggsize_ = a; };
			int particles() { return (this->atoms_ - this->cornrs_); }
			int& corners() { return this->cornrs_; }
			void corners(int c) { this->cornrs_ = c; }
			str& comment(){ return this->comment_; }
			void comment(cstr &str_){ this->comment_ = str_; }
			type& type(){ return this->type_; }
			void type(simReplay::type ty_){ this->type_ = ty_; }
			float& xmin(){ return this->xmin_; }
			float& ymin(){ return this->ymin_; }
			float& zmin(){ return this->xmin_; };
			void xmin(const float& x){ this->xmin_ = x; }
			void ymin(const float& y){ this->ymin_ = y; }
			void zmin(const float& z){ this->xmin_ = z; }
			float& xmax(){ return this->xmax_; }
			float& ymax(){ return this->ymax_; }
			float& zmax(){ return this->xmax_; };
			void xmax(const float& x){ this->xmax_ = x; }
			void ymax(const float& y){ this->ymax_ = y; }
			void zmax(const float& z){ this->xmax_ = z; }
			float hx(){ return (this->xmax_ - this->xmin_); }
			float hy(){ return (this->ymax_ - this->ymin_); }
			float hz(){ return (this->zmax_ - this->zmin_); }
			void print();
		};

		//.. print state of simreplay file properties
		void properties::print(){
			printf("\nset:\t\t%I", this->set_);
			printf("\natoms:\t\t%i", this->atoms_);
			printf("\naggsize:\t%i", this->aggsize_);
			printf("\ncorners:\t%i", this->cornrs_);
			printf("\ncomment:\t%s", this->comment_.c_str());
			printf("\ntype:\t\t%i", (int)this->type_);
			printf("\nx-range:\t%f - %f", this->xmin_, this->xmax_);
			printf("\ny-range:\t%f - %f", this->ymin_, this->ymax_);
			printf("\nz-range:\t%f - %f", this->zmin_, this->zmax_);
		}

		//.. get properties of simreplay file (not currently opened).
		//   returns false if does not complete.
		//   assumes 4 corners because of default values in properties
		bool getProperties(str& file, properties& new_props){
			std::ifstream fin(file.c_str(), std::ios::in); // open file
			str line;
			if (!fin.is_open()) return false;
			if (!std::getline(fin, line)) return false; // atom number
			new_props.atoms((int)std::strtod(line.c_str(), NULL));
			if (!std::getline(fin, line)) return false; // comment line
			new_props.comment(line);

			if (!std::getline(fin, line)) return false; // get first line
			int cols = 0; char ty_1st;
			sstrm ss(line); str tmp; 
			if (ss >> ty_1st) cols++; // grab first particle type
			while (ss >> tmp) cols++; // count extra cols
			char ty_; int aggregate = 0; 
			do {
				aggregate++; // add 1 to aggsize
				ss.clear();
				if (!getline(fin, line)) return false;
				ss.str(line);
				ss >> ty_;
			} while (ty_ == ty_1st); // until new particle type read
			new_props.aggsize(aggregate); // assign aggregate size

			if (cols == 4) new_props.type(type::xyz);
			else if (cols == 5) new_props.type(type::xyzc);
			else if (cols == 7) new_props.type(type::xyzv);
			else new_props.type(type::no_type);

			const int more_to_read = (new_props.atoms() - new_props.corners() - aggregate - 1);
			for (int i = 0; i < more_to_read; i++)
				if (!std::getline(fin, line)) return false; // loop to end
			for (int i = 0; i < new_props.corners(); i++){ // get box corners
				if (!std::getline(fin, line)) return false;
				sstrm css(line); char c; float x, y, z;
				css >> c >> x >> y >> z; // seek for box size
				new_props.xmin() = std::min(new_props.xmin(), x);
				new_props.xmax() = std::max(new_props.xmax(), x);
				new_props.ymin() = std::min(new_props.ymin(), y);
				new_props.ymax() = std::max(new_props.ymax(), y);
				new_props.zmin() = std::min(new_props.zmin(), z);
				new_props.zmax() = std::max(new_props.zmax(), z);
			}
			new_props.set();
			fin.close(); // close file
			return true; // function finished correctly
		}

		//.. read file from open ifstream
		template<typename T = float>
		bool readParticles(std::ifstream& fin, properties& props,
			T* x = NULL, T* y = NULL, T* z = NULL,
			T* c_or_vx = NULL, T* vy = NULL, T* vz = NULL){
			if (!props.is_set()) return false; // properties not setup, stop
			if (!fin.is_open()) return false; // file stream not open, stop
			
			str line;
			if (!std::getline(fin, line)) return false; // atom number
			const int atoms = (int)std::strtod(line.c_str(), NULL);
			if (atoms != props.atoms()){ // check properties agree
				printf("\nWarning: adjusting properties for frame");
				props.atoms() = atoms;
			}
			if (!std::getline(fin, line)) return false; // comment line
			if (std::strcmp(line.c_str(), props.comment().c_str()) != 0){ // adjust comment if changed
				props.comment(line);
			}

			char charTrash; float junk;
			switch (props.type()) { // read lines in proper way
			case type::xyz: // -------------------------------
				for (int i = 0; i < props.particles(); i++)
					if (!readLine(fin, charTrash, x[i], y[i], z[i]))
						return false;
				for (int i = 0; i < props.corners(); i++)
					if (!readLine(fin, charTrash, junk, junk, junk))
						return false;
				break;
			case type::xyzc: // ------------------------------------
				for (int i = 0; i < props.particles(); i++)
					if (!readLine(fin, charTrash, x[i], y[i], z[i], c_or_vx[i])) 
						return false;
				for (int i = 0; i < props.corners(); i++)
					if (!readLine(fin, charTrash, junk, junk, junk, junk))
						return false;
				break;
			case type::xyzv: // ---------------------------------------
				for (int i = 0; i < props.particles(); i++)
					if (!readLine(fin, charTrash, x[i], y[i], z[i], c_or_vx[i], vy[i], vz[i])) 
						return false;
				for (int i = 0; i < props.corners(); i++)
					if (!readLine(fin, charTrash, junk, junk, junk, junk, junk, junk))
						return false;
			case type::no_type: // -------------------------------------
				return false;
			}
			return true;
		}

		//.. implement read lines
		bool readLine(std::ifstream& fin, char& c, float& x, float& y, float& z){
			str line;
			if (!std::getline(fin, line)) return false;
			sstrm row(line);
			row >> c >> x >> y >> z;
			return true;
		}
		bool readLine(std::ifstream& fin, char& c, float& x, float& y, float& z, float& C){
			str line;
			if (!std::getline(fin, line)) return false;
			sstrm row(line);
			row >> c >> x >> y >> z >> C;
			return true;
		}
		bool readLine(std::ifstream& fin, char& c, float& x, float& y, float& z, float& vx, float& vy, float& vz){
			str line;
			if (!std::getline(fin, line)) return false;
			sstrm row(line);
			row >> c >> x >> y >> z >> vx >> vy >> vz;
			return true;
		}
	}

	//.. utilies for nematic worms only
	namespace nw {

		//.. calculate distances in periodic box
		void pbc(float &dx, const float& Lx) {
			if (dx > Lx / 2) dx -= Lx;
			else if (dx < -Lx / 2) dx += Lx;
		}

		//.. return true position in period box
		void pbc_pos(float &x, const float& Lx) {
			if (x > Lx) x -= Lx;
			else if (x < 0) x += Lx;
		}
	}
}