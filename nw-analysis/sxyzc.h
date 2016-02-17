// Nematic Worm Analysis
// 2.17.16
// Produce xyzc of order parameter
// Mike Varga
#ifndef __SXYZC_H__
#define __SXYZC_H__
// ---------------------------
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <stdio.h>
#include <math.h>
#include <vector>
#include "utilities.h"
#include "liquid_crystals.h"
// --------------------------
namespace sxyzc {

	const std::string funcName = "sxyzc";

	static void show_usage(std::string name)
	{
		std::cerr << "\a";
		std::cerr << "Usage: " << name << " <option(s)> file...\n"
			<< "\t-i,--input\t\tInput file name\n"
			<< "\t-o,--output\t\tOutput file name\n"
			<< "Options:\n"
			<< "\t-h,--help\t\tShow this help message\n"
			<< "If no in/out options specified, default output file name is 'Q-tensor.txt'\n"
			<< "and input file must follow program instance\n"
			<< std::endl;
	}

	static int process_arg(std::ifstream &fin, std::ofstream &fout,
		std::vector<std::string> argv)
	{
		const int argc = argv.size();
		for (int i = 0; i < argc; ++i){
			std::string arg = argv[i];
			if ((arg == "-h") || (arg == "--help")){
				show_usage(funcName);
				return 1;
			}
			else if ((arg == "-i") || (arg == "--input")){
				if (i + 1 < argc){
					fin.open(argv[i + 1]);
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 2;
				}
			}
			else if ((arg == "-o") || (arg == "--output")){
				if (i + 1 < argc){
					fout.open(argv[i + 1]);
					i++;
				}
				else{
					std::cerr << "--output option requires one argument." << std::endl;
					return 3;
				}
			}
			else{
				fin.open(argv[i]);
				fout.open(funcName + ".xyzc");
				return 0;
			}
		}
		return 0;
	}

	int calculate(std::vector<std::string> argv){

		std::ifstream fin;
		std::ofstream fout;

		//.. process command line sub arguments
		int process_arg_status = process_arg(fin, fout, argv);
		if (process_arg_status != 0) return process_arg_status;

		if (!fin.is_open()){
			std::cerr << "Error opening file" << std::endl
				<< "Check for correct input name and directory\n" << std::endl;
			show_usage(funcName);
			return 20;
		}
		else{
			int		numParticles;
			int		numPerWorm;
			int		numWorms;
			double	numWorms2;
			int     numFrame = 1;
			char	charTrash;
			float	k2spring;
			float	hx;
			float	hy;
			float	floatTrash;
			float	boxSize = 2.5f;

			while (!fin.eof()){
				numParticles = 4; // Deals with case of black line at end of file.
				fin >> numParticles;
				fin >> numPerWorm >> k2spring >> hx >> hy;
				numParticles -= 4;

				if (numParticles != 0){

					numWorms = numParticles / numPerWorm;
					numWorms2 = numWorms*numWorms;

					// Set up dynamic arrays 4 particles
					float** x = new float*[numWorms];
					float** y = new float*[numWorms];
					for (int w = 0; w < numWorms; w++){
						x[w] = new float[numPerWorm];
						y[w] = new float[numPerWorm];
					}

					// Determine lattice for calculating average S
					const int numXBox = int(std::ceil(hx / boxSize));
					const int numYBox = int(std::ceil(hy / boxSize));
					const int numBox = numXBox * numYBox;

					// set up dynamic array for scaler order
					float** S = new float*[numXBox];
					int**	ns = new int*[numXBox];
					for (int i = 0; i < numXBox; i++){
						S[i] = new float[numYBox];
						ns[i] = new int[numYBox];
						for (int j = 0; j < numYBox; j++){
							S[i][j] = 0;
							ns[i][j] = 0;
						}
					}

					// Read in X,Y,Z positions for frame
					std::cout << "Reading Frame " << numFrame << "\t";
					for (int w = 0; w < numWorms; w++)
						for (int p = 0; p < numPerWorm; p++)
							fin >> charTrash >> x[w][p] >> y[w][p] >> floatTrash;

					// Dump corner particles at end of file to trash
					for (int t = 0; t < 4; t++)
						fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;

					// Load S[][] with y/x for all particle pairs
					for (int w = 0; w < numWorms; w++){
						for (int p = 0; p < numPerWorm - 1; p++){
							const float Xcom = (x[w][p] + x[w][p + 1]) / 2.0f;
							const float Ycom = (y[w][p] + y[w][p + 1]) / 2.0f;
							const int	ibox = int(floor(Xcom / boxSize));
							const int	jbox = int(floor(Ycom / boxSize));

							double YdivX = (y[w][p] - y[w][p + 1]) / (x[w][p] - x[w][p + 1]);
							S[ibox][jbox] += YdivX;
							ns[ibox][jbox]++;
						}
					}

					// Get Average Y/X in box
					for (int i = 0; i < numXBox; i++){
						for (int j = 0; j < numYBox; j++){
							if (ns[i][j] != 0) S[i][j] /= float(ns[i][j]);
							else ns[i][j] = -1; // Flag empty boxes							
						}
					}

					// Handle Empty boxes
					for (int i = 0; i < numXBox; i++){
						for (int j = 0; j < numYBox; j++){
							if (ns[i][j] == -1){
								int ip1 = (i + 1) % numXBox;
								int jp1 = (j + 1) % numYBox;
								int im1 = i - 1;
								int jm1 = j - 1;
								if (i == 0) im1 = numXBox - 1;
								if (j == 0) jm1 = numYBox - 1;

								S[i][j] = (S[ip1][j] + S[i][jp1] + S[im1][j] + S[i][jm1]) / 4;
							}
						}
					}

					// Convert average Y/X to angles
					std::cout << "Analysis Completed" << std::endl;
					for (int i = 0; i < numXBox; i++)
						for (int j = 0; j < numYBox; j++)
							S[i][j] = std::atan(S[i][j]);

					fout << numBox << std::endl << "comment" << std::endl;
					for (int i = 0; i < numXBox; i++)
						for (int j = 0; j < numYBox; j++)
							fout << "A " << i << " " << j << " " << 0 << " " << std::sin(2 * S[i][j])*std::sin(2 * S[i][j]) << std::endl;

					// Handle Dynamic Arrays
					for (int dW = 0; dW < numWorms; dW++){
						delete[] x[dW];
						delete[] y[dW];
					}

					for (int di = 0; di < numXBox; di++){
						delete[] S[di];
						delete[] ns[di];
					}

					delete[] x;
					delete[] y;
					delete[] S;
					delete[] ns;
					numFrame++;
				}
			}
		}
		return EXIT_SUCCESS;
	}
//*************************************************************************

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
		int				&sfid,
		int				&efid)
	{
		const int argc = argv.size();
		for (int i = 0; i < argc; ++i)
		{
			std::string arg = argv[i];
			if ((arg == "-h") || (arg == "--help"))
			{
				show_usage(funcName);
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
			else
			{
				outbasename = funcName;
				return 0;
			}
		}
		return 0;
	}

	// -------------------------------------------------
	// Run calculation in 3d
	int calculate_3d(std::vector<std::string> argv){

		std::string finBaseName;
		std::string foutBaseName;
		float * x = { 0 };
		float * y = { 0 };
		float * z = { 0 };
		float	boxWidth = 2.0f;
		int		startFileId = -1;
		int		endFileId = -1;
		int		numFrame = 0;

		//.. process and assign cmdline args
		int process_arg_status = process_arg_3d(finBaseName,
			foutBaseName,
			argv,
			boxWidth,
			startFileId,
			endFileId);
		if (process_arg_status != 0)
			return process_arg_status;

		//.. open output files
		std::ofstream fxyz(foutBaseName + ".xyzc", std::ios::out);
		if (!fxyz.is_open()) return 10;

		//.. loop through all files (or just once when -1 & -1)
		for (int fid = startFileId; fid <= endFileId; fid++)
		{
			//.. construct file name
			std::ostringstream ifname;
			ifname << finBaseName;
			if (startFileId >= 0) ifname << fid;
			ifname << ".xyz";

			//.. gather properties and assure xyz format
			util::simReplay::properties fileProps;
			fileProps.getFrom(ifname.str());
			fileProps.print();
			if (fileProps.type != util::simReplay::type::xyz)
				return 20;

			//.. open file
			std::ifstream fin(ifname.str(), std::ios::in);
			if (!fin.is_open()){
				std::cerr << std::endl << funcName << ": Error opening file\n"
					<< funcName << ": Check for correct input name and directory";
				show_usage_3d(funcName);
				return 30;
			}

			//.. stuff
			const int numParticles = fileProps.particles();
			const int numPerWorm = fileProps.aggsize();
			const int numWorms = numParticles / numPerWorm;
			const float hx = fileProps.hx();
			const float hy = fileProps.hy();

			//.. each frame loop
			while (!fin.eof())
			{
				// allocate memory
				x = new float[numParticles];
				y = new float[numParticles];
				z = new float[numParticles];
				util::simReplay::readParticles(fin, fileProps, x, y, z);
				printf("\n%s: Frame %i read from %s", funcName.c_str(), numFrame, ifname.str().c_str());

				// container for all particle directors
				const int xdim = (int)ceil(hx / boxWidth);
				const int ydim = (int)ceil(hy / boxWidth);
				std::vector<ngroup> all_n;
				all_n.resize(xdim*ydim);

				// put directors into boxes
				float dx, dy, dz, px, py, pz;
				int i, j, b, id;
				for (int w = 0; w < numWorms; w++){
					for (int p = 0; p < numPerWorm - 1; p++){
						id = w*numPerWorm + p;
						dx = x[id + 1] - x[id];
						dy = y[id + 1] - y[id];
						dz = z[id + 1] - z[id];
						px = x[id] + dx / 2; // midpoint of vector
						py = y[id] + dy / 2;
						pz = z[id] + dz / 2;
						util::nw::pbc(dx, hx); // pbc for distances
						util::nw::pbc(dy, hy);
						util::nw::pbc_pos(px, hx); // pbc for location
						util::nw::pbc_pos(py, hy);
						i = int(px / boxWidth); // box i coord
						j = int(py / boxWidth); // box j coord
						b = j*xdim + i; // box address
						n dir(dx, dy, dz);
						all_n.at(b).push_back(dir); // add dir to box
					}
				}			

				// Calculate S and print to files
				fxyz << xdim * ydim << std::endl;
				fxyz << "Frame " << numFrame++ << std::endl;
				for (int j = 0; j < ydim; j++){
					for (int i = 0; i < xdim; i++){
						int id = j*xdim + i;
						fxyz << "A " << i << " " << j << " 0 " << calculate_S(all_n[id]) << std::endl;
					}
				}

				delete[] x;
				delete[] y;
				delete[] z;
				printf("\n%s: Frame memory deleted.", funcName.c_str());
			}
			fin.close();
			printf("\n%s: Input file closed.", funcName.c_str());
		}
		fxyz.close();
		printf("\n%s: Output files closed.", funcName.c_str());
		return EXIT_SUCCESS;
	}
}

#endif