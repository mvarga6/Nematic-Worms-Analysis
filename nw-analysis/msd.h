// Nematic Worms Analysis
// 2.4.16
// Detect thickness of worm layer
// Mike Varga
#ifndef __MSD_H__
#define __MDS_H__
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
namespace msd {
	const std::string funcName = "mean squared displacement";

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
	// Run calculation in 3d
	int calculate_3d(std::vector<std::string> argv){

		//.. calculation long variables
		std::string finBaseName;
		std::string foutBaseName;
		float * x		= { 0 };
		float * y		= { 0 };
		float * z		= { 0 };
		float * x_old	= { 0 };
		float * y_old	= { 0 };
		float * z_old	= { 0 };
		float * X		= { 0 };
		float * Y		= { 0 };
		float * Z		= { 0 };
		//float * X_old	= { 0 };
		//float * Y_old	= { 0 };
		//float * Z_old	= { 0 };
		bool set = false;
		int		startFileId = -1;
		int		endFileId = -1;
		int		numFrame = 0;

		//.. process and assign cmdline args
		int process_arg_status = process_arg_3d(finBaseName,
			foutBaseName,
			argv,
			startFileId,
			endFileId);
		if (process_arg_status != 0) return process_arg_status;

		std::vector< std::vector<float> > X0;
		std::vector< std::vector<float> > Y0;
		std::vector< std::vector<float> > Z0;

		//.. open output file
		std::ofstream fcsv(foutBaseName + ".csv", std::ios::out);
		if (!fcsv.is_open()) return 10;

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

			//.. stuff from file
			const int numParticles = fileProps.particles();
			const float hx = fileProps.hx();
			const float hy = fileProps.hy();

			if (!set){ // alloc on first time only
				X = new float[numParticles];
				Y = new float[numParticles];
				Z = new float[numParticles];
				/*X_old = new float[numParticles];
				Y_old = new float[numParticles];
				Z_old = new float[numParticles];*/
				x = new float[numParticles];
				y = new float[numParticles];
				z = new float[numParticles];
				x_old = new float[numParticles];
				y_old = new float[numParticles];
				z_old = new float[numParticles];
				set = true;
			}

			//.. each frame loop
			while (!fin.eof())
			{
				// get positions
				util::simReplay::readParticles(fin, fileProps, x, y, z);
				printf("\n%s: Frame %i read from %s", funcName.c_str(), ++numFrame, ifname.str().c_str());
				
				if (numFrame == 1){ // init to inital positions
					for (int i = 0; i < numParticles; i++){
						X[i] = x[i];
						Y[i] = y[i];
						Z[i] = z[i];
						x_old[i] = x[i];
						y_old[i] = y[i];
						z_old[i] = z[i];
					}
				}
				else{ // add local displacement to absolute positions
					float dx, dy, dz;
					for (int i = 0; i < numParticles; i++){
						dx = x[i] - x_old[i]; // local displacement
						dy = y[i] - y_old[i];
						dz = z[i] - z_old[i];
						util::nw::pbc(dx, hx);
						util::nw::pbc(dy, hy);
						X[i] += dx; // add to absolute positions
						Y[i] += dy;
						Z[i] += dz;
						x_old[i] = x[i]; // save local positions
						y_old[i] = y[i];
						z_old[i] = z[i];
					}
				}
				
				// push back new absolute position
				std::vector<float> add_X0;
				std::vector<float> add_Y0;
				std::vector<float> add_Z0;
				add_X0.reserve(numParticles);
				add_Y0.reserve(numParticles);
				add_Z0.reserve(numParticles);
				for (int i = 0; i < numParticles; i++){
					add_X0.push_back(X[i]);
					add_Y0.push_back(Y[i]);
					add_Z0.push_back(Z[i]);
				}
				X0.push_back(add_X0);
				Y0.push_back(add_Y0);
				Z0.push_back(add_Z0);
				add_X0.clear();
				add_Y0.clear();
				add_Z0.clear();

				// calculate msd
				float dX, dY, dZ;
				std::vector<float> msd;
				for (int t0 = X0.size() - 1; t0 >= 0; t0--){ // loop backwards in time
					float R2t0 = 0;
					for (int i = 0; i < numParticles; i++){
						dX = X[i] - X0[t0][i]; // distance to pos at t0 from current
						dY = Y[i] - Y0[t0][i];
						dZ = Z[i] - Z0[t0][i];
						//util::nw::pbc(dx, hx);
						//util::nw::pbc(dy, hy);
						R2t0 += dX*dX + dY*dY + dZ*dZ;
					}
					R2t0 /= (float)numParticles; // <(r(t) - r(t0))^2> ensemble ave
					msd.push_back(R2t0);
				}

				// Print current data to file
				fcsv << numFrame;
				for (auto it : msd) fcsv << ", " << it;
				fcsv << std::endl;
			}
			fin.close();
			printf("\n%s: Input file closed.", funcName.c_str());
			
		}
		delete[] x;
		delete[] y;
		delete[] z;
		delete[] X;
		delete[] Y;
		delete[] Z;
		/*delete[] X_old;
		delete[] Y_old;
		delete[] Z_old;*/
		delete[] x_old;
		delete[] y_old;
		delete[] z_old;
		fcsv.close();
		printf("\n%s: Output files closed.", funcName.c_str());
		return EXIT_SUCCESS;
	}
}

#endif