#pragma once
// Nematic Worm Analysis
// 2.9.16
// Spatial Correlation Functions 2D and 3D
// Mike Varga
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "utilities.h"

namespace gofr {

	const float _PI = 3.14159265359f;
	const float _2PI = 2 * _PI;

	const std::string funcName = "gofr";
	//------------------------------------------------------
	static void show_usage_2d(std::string name)
	{
		std::cerr << "\a";
		std::cerr << "Usage: " << name << " <option(s)> file...\n"
			<< "\t-i,--input\t\tInput file name\n"
			<< "\t-o,--output\t\tOutput file name\n"
			<< "Options:\n"
			<< "\t-h,--help\t\tShow this help message\n"
			<< "\t-m,--max\t\tMaximum correlation distance (default=10.0)\n"
			<< "\t-bw,--width\t\tMaximum correlation distance (default=10.0)\n"
			<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
			<< "\t-e,--end\t\tDefine id of last input file (default=1)\n"
			<< "If no in/out options specified, default output file name is 'Q-tensor.txt'\n"
			<< "and input file must follow program instance\n"
			<< std::endl;
	}
	
	//------------------------------------------------------
	static int process_arg_2d(std::string &basename,
		std::ofstream 	&fout,
		std::vector<std::string> argv,
		float 			&range,
		float			&binwidth,
		int				&sfid,
		int				&efid)
	{
		const int argc = argv.size();
		for (int i = 0; i < argc; ++i){
			std::string arg = argv[i];
			if ((arg == "-h") || (arg == "--help")){
				show_usage_2d(funcName);
				return 1;
			}
			else if ((arg == "-i") || (arg == "--input")){
				if (i + 1 < argc){
					basename = argv[i + 1];
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 2;
				}
			}
			else if ((arg == "-m") || (arg == "--max")){
				if (i + 1 < argc){
					range = strtof(argv[i + 1].c_str(), NULL);
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 3;
				}
			}
			else if ((arg == "-bw") || (arg == "--width")){
				if (i + 1 < argc){
					binwidth = strtof(argv[i + 1].c_str(), NULL);
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 4;
				}
			}
			else if ((arg == "-s") || (arg == "--start")){
				if (i + 1 < argc){
					sfid = int(strtod(argv[i + 1].c_str(), NULL));
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 5;
				}
			}
			else if ((arg == "-e") || (arg == "--end")){
				if (i + 1 < argc){
					efid = int(strtod(argv[i + 1].c_str(), NULL));
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 6;
				}
			}
			else if ((arg == "-o") || (arg == "--output")){
				if (i + 1 < argc){
					fout.open(argv[i + 1]);
					i++;
				}
				else{
					std::cerr << "--output option requires one argument." << std::endl;
					return 7;
				}
			}
			else{
				fout.open(funcName + ".csv");
				return 0;
			}
		}
		return 0;
	}

	/*************************************************
	*	Two-dimensional (xy-plane) g(r) calculation 
	***************************************************/
	int calculate_2d(std::vector<std::string> argv){

		//.. i/o things
		std::string finBaseName;
		std::ofstream fout;

		//.. user input parameters
		float intRange = 10.0f;
		int	   startFileId = -1;
		int	   endFileId = -1;
		int     numFrame = 1;
		float	binWidth = 1.0f;

		//.. process command line sub arguments
		int process_arg_status = process_arg_2d(finBaseName, fout, argv, intRange, binWidth, startFileId, endFileId);
		if (process_arg_status != 0) return process_arg_status;

		//.. number of bins and first line of data file
		int numBin = int(ceil(intRange / binWidth));
		fout << "bins";
		for (int i = 1; i <= numBin; i++)
			fout << ", " << float(i)*binWidth;
		fout << std::endl;

		//.. linked list parameters (for worms COM, not particles)
		float dcell = float(intRange);
		int nxcell, nycell, ncell;
		const int ddx[5] = { 0, -1, 0, 1, 1 };
		const int ddy[5] = { 0, 1, 1, 1, 0 };

		//.. loop through all input files
		for (int fid = startFileId; fid <= endFileId; fid++)
		{
			//.. make ifstream and open
			std::ostringstream ifname;
			ifname << finBaseName;
			if (startFileId >= 0) ifname << fid;
			ifname << ".xyz";
			std::ifstream fin(ifname.str(), std::ios::in);
			if (!fin.is_open()){
				std::cerr << std::endl << funcName << ": Error opening file\n"
					<< funcName << ": Check for correct input name and directory";
				show_usage_2d(funcName);
				return 20;
			}

			// Simulation parameters
			int		numParticles;
			int		numPerWorm;
			int		numWorms;
			float	numWorms2;
			char	charTrash;
			float	k2spring;
			float	hx;
			float	hy;
			float	floatTrash;
			float	intRange2 = intRange * intRange;

			//.. loop through entire file
			while (!fin.eof())
			{
				numParticles = 4; // Deals with case of black line at end of file.
				fin >> numParticles;
				fin >> numPerWorm >> k2spring >> hx >> hy;
				numParticles -= 4;

				if (numParticles == 0) break;

				numWorms = numParticles / numPerWorm;
				numWorms2 = numWorms*numWorms;
				const float hxo2 = hx / 2.0;
				const float hyo2 = hy / 2.0;

				//.. make linked list
				nxcell = int(ceil(hx / dcell));
				nycell = int(ceil(hy / dcell));
				ncell = nxcell*nycell;

				// Allocate Containers
				float * x = new float[numParticles];
				float * y = new float[numParticles];
				int * ptz = new int[numParticles];
				int * heads = new int[ncell];

				//.. init heads
				for (int i = 0; i < ncell; i++)
					heads[i] = -1;

				// Read in X,Y,Z positions for frame and set linked list 
				std::cout << "Reading frame " << numFrame << " from " << ifname.str() << std::endl;
				for (int i = 0; i < numParticles; i++)
				{
					fin >> charTrash >> x[i] >> y[i] >> floatTrash;

					//.. put into linked list
					int icell = int(x[i] / dcell);
					int jcell = int(y[i] / dcell);
					int scell = jcell*nxcell + icell;
					ptz[i] = heads[scell];
					heads[scell] = i;
				}

				// Dump corner particles at end of file to trash
				for (int t = 0; t < 4; t++)
					fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;

				// Components of calculation
				std::vector<float> G;
				float aveOver = 1;
				for (int i = 0; i <= numBin; i++)
					G.push_back(0.0f);

				//.. calculate using linked list
				std::cout << "Calculating ... \n";
				for (int ic = 0; ic < nxcell; ic++)
				{
					for (int jc = 0; jc < nycell; jc++)
					{
						//.. scalar add of cell
						int scell = jc*nxcell + ic;

						//.. look for stop flag
						if (heads[scell] == -1) continue;

						//.. loop over adjacent cells
						for (int dir = 0; dir < 5; dir++)
						{
							//.. neighbor cell scalor address
							int icnab = (ic + ddx[dir]) % nxcell;
							int jcnab = (jc + ddy[dir]) % nycell;
							if (icnab < 0) icnab += nxcell;
							if (jcnab < 0) jcnab += nycell;
							int scnab = jcnab * nxcell + icnab;

							//.. loop for stop flag
							if (heads[scnab] == -1) continue;

							int ii = heads[scell];
							while (ii >= 0)
							{
								int jj = heads[scnab];
								while (jj >= 0)
								{
									//.. calculate distance
									float dx = x[ii] - x[jj];
									float dy = y[ii] - y[jj];
									if (dx > hxo2) dx -= hx;
									if (dx < -hxo2) dx += hx;
									if (dy > hyo2) dy -= hy;
									if (dy < -hyo2) dy += hy;
									float r2 = dx*dx + dy*dy;

									//.. MAIN CALCULATION
									if (r2 <= intRange2)
									{
										//.. find bin
										float r = sqrt(r2);
										int b = int(r / binWidth);
										if ((b >= 0) && (b <= numBin))
										{
											//.. count
											aveOver += 1;
											G.at(b) += 1;
										}
									}
									jj = ptz[jj];
								}
								ii = ptz[ii];
							}
						}
					}
				}

				//.. Average G(r) and print to out file
				fout << numFrame++;
				for (int i = 0; i < numBin; i++)
				{
					float r1 = float(i)*binWidth;
					float r2 = float(i + 1)*binWidth;
					float area = _PI*(r2*r2 - r1*r1);
					G[i] /= aveOver;
					G[i] *= 2.0f;
					G[i] /= area;
					fout << ", " << G[i];
				}
				fout << std::endl;

				//.. eliminate dynamic memory
				delete[] x;
				delete[] y;
				delete[] ptz;
				delete[] heads;
			} // End of File

			fin.close();
		}
		fout.close();
		return EXIT_SUCCESS;
	}
//*******************************************************************************************
	static void show_usage_3d(std::string name)
	{
		std::cerr << "\a";
		std::cerr << "Usage: " << name << " <option(s)> file...\n"
			<< "\t-i,--input\t\tInput file name\n"
			<< "\t-o,--output\t\tOutput file name\n"
			<< "Options:\n"
			<< "\t-h,--help\t\tShow this help message\n"
			<< "\t-m,--max\t\tMaximum correlation distance (default=10.0)\n"
			<< "\t-bw,--width\t\tG(r) resolution (default=1.0)\n"
			<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
			<< "\t-e,--end\t\tDefine id of last input file (default=1)\n"
			<< "If no in/out options specified, default output file name is 'Q-tensor.txt'\n"
			<< "and input file must follow program instance\n"
			<< std::endl;
	}
	
	// ---------------------------------------------------
	static int process_arg_3d(std::string &basename,
		std::string 	&outbasename,
		std::vector<std::string> argv,
		float 			&range,
		float			&binwidth,
		int				&sfid,
		int				&efid)
	{
		const int argc = argv.size();
		for (int i = 0; i < argc; ++i){
			std::string arg = argv[i];
			if ((arg == "-h") || (arg == "--help")){
				show_usage_3d(funcName);
				return 1;
			} // --------------------------------------------
			else if ((arg == "-i") || (arg == "--input")){
				if (i + 1 < argc){
					basename = argv[i + 1];
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 2;
				}
			} // --------------------------------------------
			else if ((arg == "-o") || (arg == "--output")){
				if (i + 1 < argc){
					outbasename = argv[i + 1];
					i++;
				}
				else{
					std::cerr << "--output option requires one argument." << std::endl;
					return 3;
				}
			} // --------------------------------------------
			else if ((arg == "-m") || (arg == "--max")){
				if (i + 1 < argc){
					range = strtof(argv[i + 1].c_str(), NULL);
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 4;
				}
			} // --------------------------------------------
			else if ((arg == "-bw") || (arg == "--width")){
				if (i + 1 < argc){
					binwidth = strtof(argv[i + 1].c_str(), NULL);
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 5;
				}
			} // --------------------------------------------
			else if ((arg == "-s") || (arg == "--start")){
				if (i + 1 < argc){
					sfid = int(strtod(argv[i + 1].c_str(), NULL));
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 6;
				}
			} // --------------------------------------------
			else if ((arg == "-e") || (arg == "--end")){
				if (i + 1 < argc){
					efid = int(strtod(argv[i + 1].c_str(), NULL));
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 7;
				}
			} // --------------------------------------------
			else{
				outbasename = funcName;
				return 0;
			}
		}
		return 0;
	}

	/*************************************************
	*	Three-dimensional (xyz-space) g(r) calculation
	***************************************************/
	int calculate_3d(std::vector<std::string> argv){

		std::string finBaseName; //input file base name
		std::string foutBaseName(funcName); // output file base name
		float * x = { 0 }; // x positions
		float * y = { 0 }; // y positions
		float * z = { 0 }; // z positions
		int * ptz = { 0 }; // linked list pointer
		int * heads = { 0 }; // head of each cell
		const int ddx[5] = { 0, -1, 0, 1, 1 }; // vector comps to next box
		const int ddy[5] = { 0, 1, 1, 1, 0 };
		float	binWidth = 1.0f; // g(r) bin size
		float	range = 20.0f; // maximum distance
		int		startFileId = -1; // start file id
		int		endFileId = -1; // end file id
		int		numFrame = 1; // frame count
		
		std::vector<float> G; // correlation functions
		long long int sum = 0;
		int numBin = int(range / binWidth);
		for (int i = 0; i <= numBin; i++)
			G.push_back(0.0f); // init to zero
		numBin = G.size();

		//.. process and assign cmdline args
		int process_arg_status = process_arg_3d(finBaseName,
			foutBaseName,
			argv,
			range,
			binWidth,
			startFileId,
			endFileId);
		if (process_arg_status != 0)
			return process_arg_status;

		// open output and over-write
		std::ofstream fcsv(foutBaseName + ".csv", std::ios::out);
		if (!fcsv.is_open()) // check for file being open
			return 10;

		//.. loop through all files (or just once when -1 & -1)
		for (int fid = startFileId; fid <= endFileId; fid++)
		{
			//.. construct name
			std::ostringstream ifname;
			ifname << finBaseName; 
			if (startFileId >= 0) ifname << fid; // add number if needed
			ifname << ".xyz";
			util::simReplay::properties fileProps; // file proerties
			fileProps.getFrom(ifname.str()); // gather properties
			fileProps.print(); // show properties

			//.. open file for reading
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
			const float hxo2 = hx / 2.0f;
			const float hyo2 = hy / 2.0f;
			const int nxcell = int(ceil(hx / range));
			const int nycell = int(ceil(hy / range));
			const int ncell = nxcell*nycell;

			// allocate memory
			x = new float[numParticles];
			y = new float[numParticles];
			z = new float[numParticles];
			ptz = new int[numParticles];
			heads = new int[ncell];
			printf("\nMemory allocated.");

			//.. each frame loop
			//while (!fin.eof())
			while (util::simReplay::readParticles(fin, fileProps, x, y, z))
			{
				//. reset neighbors list
				for (int i = 0; i < ncell; i++) 
					heads[i] = -1; //.. init heads

				//if (!util::simReplay::readParticles(fin, fileProps, x, y, z)) break;
				printf("\n%s: Frame %i read from %s", funcName.c_str(), numFrame++, ifname.str().c_str());

				//.. construct linked list
				for (int i = 0; i < numParticles; i++){
					int icell = int(x[i] / range);
					int jcell = int(y[i] / range);
					int scell = jcell*nxcell + icell;
					ptz[i] = heads[scell];
					heads[scell] = i;
				}

				//.. puts counts in g(r) using linked list
				printf("\nCalculating ...\t\t\t");
				float dx, dy, dz, r, rr;
				sum = 0;
				const float range2 = range*range;
				for (int ic = 0; ic < nxcell; ic++){
					for (int jc = 0; jc < nycell; jc++){
						int scell = jc*nxcell + ic; // scalar add of cell
						if (heads[scell] == -1) continue;// look for stop flag
						for (int dir = 0; dir < 5; dir++) { // loop over adjacent cells
							int icnab = (ic + ddx[dir]) % nxcell; // neighbor cell scalor address
							int jcnab = (jc + ddy[dir]) % nycell;
							if (icnab < 0) icnab += nxcell;
							if (jcnab < 0) jcnab += nycell;
							int scnab = jcnab * nxcell + icnab;
							if (heads[scnab] == -1) continue; //.. loop for stop flag
							int ii = heads[scell];
							while (ii >= 0){
								int jj = heads[scnab];
								while (jj >= 0){
									dx = x[ii] - x[jj]; // displacements
									dy = y[ii] - y[jj];
									dz = z[ii] - z[jj];
									if (dx > hxo2) dx -= hx; // x pbc
									if (dx < -hxo2) dx += hx;
									if (dy > hyo2) dy -= hy; // y pbc
									if (dy < -hyo2) dy += hy;
									rr = dx*dx + dy*dy + dz*dz; // dist^2
									if (rr <= range2){ // add to g(r)
										r = sqrt(rr); // find bin
										int b = int(r / binWidth);
										if ((b < 0) || (b > numBin)) continue; 
										G.at(b) += 1.0f; sum++;
									}
									jj = ptz[jj];
								}
								ii = ptz[ii];
							}
						}
					}
				}
				printf("done");

				// open output and over-write
				//std::ofstream fcsv(foutBaseName + ".csv", std::ios::app);
				//if (!fcsv.is_open()) // check for file being open
					//return 10;

				//.. write current average to file
				printf("\nAppending output file ...\t");
				float r1, r2, area;
				const int length = G.size();
				for (int i = 0; i < length; i++){
					r1 = float(i)*binWidth;
					r2 = float(i + 1)*binWidth;
					area = _PI*(r2*r2 - r1*r1);
					G.at(i) /= (float)sum;
					G.at(i) *= 2.0f;
					G.at(i) /= area;
					G.at(i) /= G.at(0); // starts at 1
					fcsv << G.at(i); // print value to file
					if (i != length - 1) // dont print comma in last col
						fcsv << ", ";
				}
				fcsv << std::endl;
				//fcsv.close(); // close file for until next frame
				printf("done");
				
				
			}
			fin.close();
			printf("\n%s: Input file closed.", funcName.c_str());

			printf("\nClearing memory ... ");
			delete[] x; printf("x ");
			delete[] y; printf("y ");
			delete[] z; printf("z ");
			delete[] ptz; printf("ptz ");
			delete[] heads; printf("heads ");
			printf("done");
		}
		fcsv.close();
		printf("\n%s: Output files closed.", funcName.c_str());
		return EXIT_SUCCESS;
	}
}