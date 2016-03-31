// Nematic Worms Analysis
// 2.4.16
// Detect thickness of worm layer
// Mike Varga
#ifndef __VELOCITY_MAP_H__
#define __VELOCITY_MAP_H__
// -------------------------------
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include "utilities.h"
// -------------------------------
namespace velocity {
	const std::string funcName = "velcity-map";

	// --------------------------------------------
	// Prints to user the program usage.
	static void show_usage_3d(std::string name)
	{
		std::cerr << "\a";
		std::cerr << "Usage: " << name << " <option(s)> file...\n"
			<< "\t-i,--input\t\tAdd input file (w/o ext)\n"
			<< "\t-o,--output\t\tOutput file name (w/o ext)\n"
			<< "Options:\n"
			<< "\t-h,--help\t\tShow this help message\n"
			<< "\t-w,--boxwidth\t\tSet width of integration cells\n"
			<< "\t-d,--dimension\t\tSet normal of grid plane (default=2:z-axis)\n"
			<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
			<< "\t-e,--end\t\tDefine id of last input file (default=1)\n"
			<< "If no in/out options specified, default output file name is 'unnamed.txt'\n"
			<< "and input file must follow program instance\n"
			<< std::endl;
	}

	// ------------------------------------------------
	// Process cmdline args and assign as needed.
	static int process_arg_3d(std::string		&basename,
		std::string		&outbasename,
		std::vector<std::string> argv,
		float			&boxwidth,
		float			&del_t,
		bool			&pbc,
		bool			&normal,
		int				&sfid,
		int				&efid)
	{
		const int argc = argv.size();
		for (int i = 0; i < argc; ++i)
		{
			std::string arg = argv[i];
			if ((arg == "-h") || (arg == "--help"))
			{
				show_usage_3d(funcName);
				return 1;
			}
			// ------------------------------------------
			else if ((arg == "-i") || (arg == "--input"))
			{
				if (i + 1 < argc){
					basename = argv[i + 1];
					i++;
				}
				else
				{
					std::cerr << "--input option requires one argument." << std::endl;
					return 2;
				}
			}
			// ------------------------------------------
			else if ((arg == "-o") || (arg == "--output"))
			{
				if (i + 1 < argc){
					outbasename = argv[i + 1];
					i++;
				}
				else
				{
					std::cerr << "--output option requires one argument." << std::endl;
					return 3;
				}
			}
			// ------------------------------------------
			else if ((arg == "-pbc") || (arg == "--periodic"))
			{
				pbc = true;
			}
			// ------------------------------------------
			else if ((arg == "-norm") || (arg == "--normalize"))
			{
				normal = true;
			}
			// ------------------------------------------
			else if ((arg == "-w") || (arg == "--boxwidth"))
			{
				if (i + 1 < argc){
					float assign = strtof(argv[i + 1].c_str(), NULL);
					boxwidth = assign;
					i++;
				}
				else
				{
					std::cerr << "--output option requires one argument." << std::endl;
					return 4;
				}
			}
			// ------------------------------------------
			else if ((arg == "-dt") || (arg == "--time-int"))
			{
				if (i + 1 < argc){
					float assign = strtof(argv[i + 1].c_str(), NULL);
					del_t = assign;
					i++;
				}
				else
				{
					std::cerr << "--output option requires one argument." << std::endl;
					return 4;
				}
			}
			// ------------------------------------------
			else if ((arg == "-s") || (arg == "--start")){
				if (i + 1 < argc){
					sfid = int(strtod(argv[i + 1].c_str(), NULL));
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 6;
				}
			}
			// ------------------------------------------
			else if ((arg == "-e") || (arg == "--end")){
				if (i + 1 < argc){
					efid = int(strtod(argv[i + 1].c_str(), NULL));
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 7;
				}
			}
			// ------------------------------------------
			else {
				outbasename = funcName;
				return 0;
			}
		}
		return 0;
	}

	// -------------------------------------------------
	// calculate root mean square displacement from layer
	// in uv-plane (rms of Z(u,v) )
	void vel_of_layer(float *u, float *v,
		float *Z, const float& layerZ,
		const int nparts, const float& boxwidth,
		const int udim, const int vdim,
		float **result){

		//.. tmp memory for counting in each square
		int **cnts = new int*[udim];
		for (int i = 0; i < udim; i++){
			cnts[i] = new int[vdim];
			for (int j = 0; j < vdim; j++){
				cnts[i][j] = 1;
			}
		}

		//.. loop and calcuate r2 sum
		int b1, b2;
		for (int i = 0; i < nparts; i++){
			b1 = (int)(u[i] / boxwidth);
			b2 = (int)(v[i] / boxwidth);

			if (b1 < 0 || b2 < 0) continue;
			if (b1 >= udim || b2 >= vdim) continue;

			//.. calculate rms from layer
			result[b1][b2] += (Z[i] - layerZ)*(Z[i] - layerZ);
			cnts[b1][b2]++;
		}

		//.. turn result into sqrt of average in box
		for (int i = 0; i < udim; i++){
			for (int j = 0; j < vdim; j++){
				float ans = result[i][j] / float(cnts[i][j]);
				ans = std::sqrt(ans);
				result[i][j] = ans;
			}
		}

		//.. free tmp memory
		for (int i = 0; i < udim; i++){
			delete[] cnts[i];
		}
		delete[] cnts;
	}

	// -------------------------------------------------
	// Run calculation in 3d
	int calculate_3d(std::vector<std::string> argv){

		std::string finBaseName;
		std::string foutBaseName;
		float * x = { 0 };
		float * y = { 0 };
		float * z = { 0 };
		float * x_old = { 0 };
		float * y_old = { 0 };
		float * z_old = { 0 };
		float	boxWidth = 2.0f;
		float	dt = 1.0f;
		int		startFileId = -1;
		int		endFileId = -1;
		int		numFrame = 0;
		bool	pbc = false;
		bool	norm = false;

		//.. process and assign cmdline args
		int process_arg_status = process_arg_3d(finBaseName,
			foutBaseName,
			argv,
			boxWidth,
			dt,
			pbc,
			norm,
			startFileId,
			endFileId);
		if (process_arg_status != 0)
			return process_arg_status;

		//.. open output files
		std::ofstream fxyz(foutBaseName + ".xyzv", std::ios::out);
		std::ofstream fcsv(foutBaseName + ".csv", std::ios::out);
		if (!fxyz.is_open() || !fcsv.is_open())
			return 10;

		//.. loop through all files (or just once when -1 & -1)
		for (int fid = startFileId; fid <= endFileId; fid++)
		{
			//.. construct file name
			std::ostringstream ifname;
			ifname << finBaseName;
			if (startFileId >= 0) ifname << fid;
			ifname << ".xyz";

			//.. gather properties
			util::simReplay::properties fileProps;
			fileProps.getFrom(ifname.str());
			fileProps.print();

			//.. open file
			std::ifstream fin(ifname.str(), std::ios::in);
			if (!fin.is_open()){
				std::cerr << std::endl << funcName << ": Error opening file\n"
					<< funcName << ": Check for correct input name and directory";
				show_usage_3d(funcName);
				return 20;
			}

			//.. stuff
			const int numParticles = fileProps.particles();
			const float hx = fileProps.hx();
			const float hy = fileProps.hy();

			// allocate memory
			x = new float[numParticles];
			y = new float[numParticles];
			z = new float[numParticles];
			x_old = new float[numParticles];
			y_old = new float[numParticles];
			z_old = new float[numParticles];

			//.. each frame loop
			while (!fin.eof())
			{
				util::simReplay::readParticles(fin, fileProps, x, y, z);
				printf("\n%s: Frame %i read from %s", funcName.c_str(), numFrame, ifname.str().c_str());

				if (numFrame++ == 0){
					for (int i = 0; i < numParticles; i++){
						x_old[i] = x[i];
						y_old[i] = y[i];
						z_old[i] = z[i];
					}
					continue;
				}

				// create velocity grid and counts
				const int xdim = (int)ceil(hx / boxWidth);
				const int ydim = (int)ceil(hy / boxWidth);
				float *** vel_ = new float **[xdim];
				int ** cnts_ = new int *[xdim];
				for (int i = 0; i < xdim; i++){
					vel_[i] = new float*[ydim];
					cnts_[i] = new int[ydim];
					for (int j = 0; j < ydim; j++){
						vel_[i][j] = new float[3]; // vels in 3d, grid 2d
						cnts_[i][j] = 0;
						for (int k = 0; k < 3; k++)
							vel_[i][j][k] = 0.0f;
					}
				}

				// calculate and add to velocity bins
				int xbox, ybox;
				float dx, dy, dz;
				for (int i = 0; i < numParticles; i++){
					//.. grid location
					xbox = (int)(x[i] / boxWidth);
					ybox = (int)(y[i] / boxWidth);

					//.. ignore anything outside of box
					if (xbox < 0 || xbox >= xdim) continue;
					if (ybox < 0 || ybox >= ydim) continue;

					dx = (x[i] - x_old[i]);
					dy = (y[i] - y_old[i]);
					dz = (z[i] - z_old[i]);

					if (pbc) { // adjust distance for pbc if needed
						if (dx > hx / 2) dx -= hx;
						if (dx < -hx / 2) dx += hx;
						if (dy > hy / 2) dy -= hy;
						if (dy < -hy / 2) dy += hy;
					}

					// assign velocities
					vel_[xbox][ybox][0] += dx / dt;
					vel_[xbox][ybox][1] += dy / dt;
					vel_[xbox][ybox][2] += dz / dt;

					cnts_[xbox][ybox]++; // count

					//.. save to old
					x_old[i] = x[i];
					y_old[i] = y[i];
					z_old[i] = z[i];
				}

				// Print to files
				fxyz << xdim * ydim << std::endl;
				fxyz << "Frame " << numFrame << std::endl;
				float vx, vy, vz; int cnts;
				for (int j = 0; j < ydim; j++){
					for (int i = 0; i < xdim; i++){

						//.. calculate value to print
						cnts = cnts_[i][j];
						if (cnts > 0){
							vx = vel_[i][j][0] / (float)cnts;
							vy = vel_[i][j][1] / (float)cnts;
							vz = vel_[i][j][2] / (float)cnts;
							if (norm){
								float mag_ = sqrt(vx*vx + vy*vy + vz*vz);
								vx /= mag_;
								vy /= mag_;
								vz /= mag_;
							}
						}
						else { vx = vy = vz = 0.0f; }

						//.. print xyz file
						fxyz << "A " << i << " " << j << " 0 " 
							<< vx << " " << vy << " " << vz << std::endl;
						//.. print csv file
						fcsv << sqrt(vel_[i][j][0] * vel_[i][j][0] + vel_[i][j][1] * vel_[i][j][1] + vel_[i][j][2] * vel_[i][j][2]) << ", ";
					}
					fcsv << std::endl;
				}

				for (int i = 0; i < xdim; i++){
					for (int j = 0; j < ydim; j++)
						delete[] vel_[i][j];
					delete[] vel_[i];
				} delete[] vel_;
				printf("\n%s: Frame memory deleted.", funcName.c_str());
			}
			delete[] x, y, z;
			delete[] x_old, y_old, z_old;
			fin.close();
			printf("\n%s: Input file closed.", funcName.c_str());
		}
		fxyz.close();
		fcsv.close();
		printf("\n%s: Output files closed.", funcName.c_str());
		return EXIT_SUCCESS;
	}
}

#endif