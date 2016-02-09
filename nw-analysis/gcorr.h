#pragma once
// Nematic Worm Analysis
// 11.11.14
// Theta correlation function
// Mike Varga
// ----------------------------
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
// --------------------------------------------------------------------------------
// This encapsolates a spacial correlation function calculation and things needed
// for easy use.
namespace gcorr {

	// ------------------------------
	const float _PI = 3.14159265359f;
	const float _2PI = 2 * _PI;
	// ------------------------------------
	const std::string funcName = "gcorr";
	// -----------------------------------------------
	static void show_usage(std::string name)
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
	// -----------------------------------------------
	static int process_arg(std::string &basename,
		std::ofstream 	&fout,
		std::vector<std::string> argv,
		float 			&range,
		float			&binwidth,
		int				&sfid,
		int				&efid)
	{
		int argc = argv.size();
		for (int i = 0; i < argc; ++i){
			std::string arg = argv[i];
			if ((arg == "-h") || (arg == "--help")){
				show_usage(funcName);
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

	/*************************************************************************************************
	**************************************   MAIN FUNCTION   *****************************************
	**************************************************************************************************/

	int calculate(std::vector<std::string> argv){

		//.. i/o things ************************************************************************

		std::string finBaseName;
		std::ofstream fout;

		//.. user input parameters *************************************************************

		float intRange = 10.0f;
		int	   startFileId = 1;
		int	   endFileId = 1;
		int     numFrame = 1;
		float	binWidth = 1.0f;

		//.. process command line sub arguments
		int process_arg_status = process_arg(finBaseName, 
			fout, 
			argv, 
			intRange, 
			binWidth, 
			startFileId, 
			endFileId);
		if (process_arg_status != 0){
			return process_arg_status;
		}

		//.. number of bins and first line of data file ****************************************

		int numBin = int(ceil(intRange / binWidth));
		fout << "bins";
		for (int i = 1; i <= numBin; i++)
			fout << ", " << float(i)*binWidth;
		fout << std::endl;

		//.. linked list parameters (for worms COM, not particles) *****************************

		float dcell = float(intRange);
		int nxcell, nycell, ncell;
		const int ddx[5] = { 0, -1, 0, 1, 1 };
		const int ddy[5] = { 0, 1, 1, 1, 0 };

		//.. loop through all input files ******************************************************

		for (int fid = startFileId; fid <= endFileId; fid++)
		{
			//.. make ifstream and open ********************************************************
			std::ostringstream ifname;
			ifname << finBaseName << fid << ".xyz";
			std::ifstream fin(ifname.str(), std::ios::in);
			if (!fin.is_open()){
				std::cerr << "Error opening file '" << ifname.str() << "'\n"
					<< "Check for correct input name and/or directory\n" << std::endl;
				show_usage(funcName);
				return 20;
			}

			// Simulation parameters **********************************************************

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

			//.. loop through entire file *****************************************************

			while (!fin.eof())
			{
				numParticles = 4; // Deals with case of black line at end of file.
				fin >> numParticles;
				fin >> numPerWorm >> k2spring >> hx >> hy;
				numParticles -= 4;

				if (numParticles == 0) break;

				numWorms = numParticles / numPerWorm;
				numWorms2 = numWorms*numWorms;
				float hxo2 = hx / 2.0;
				float hyo2 = hy / 2.0;

				//.. make linked list ********************************************************
				nxcell = int(ceil(hx / dcell));
				nycell = int(ceil(hy / dcell));
				ncell = nxcell*nycell;

				// Allocate Containers *******************************************************

				float * x = new float[numParticles];
				float * y = new float[numParticles];
				float * theta = new float[numParticles];
				int * ptz = new int[numParticles];
				int * heads = new int[ncell];

				//.. init heads **************************************************************

				for (int i = 0; i < ncell; i++)
					heads[i] = -1;

				// Read in X,Y,Z positions for frame and set linked list *********************
				std::cout << "Reading frame " << numFrame 
					<< " from " << ifname.str() << std::endl;
				for (int w = 0; w < numWorms; w++)
				{
					//.. read in postions for worm w *****************************************

					for (int p = 0; p < numPerWorm; p++)
					{
						int i = w*numPerWorm + p;
						fin >> charTrash >> x[i] >> y[i] >> floatTrash;
					}

					//.. calculate theta and set LL cells ************************************

					for (int p = 0; p < numPerWorm; p++)
					{
						int i = w*numPerWorm + p;
						if (p < numParticles - 1)
						{
							float dx = x[i + 1] - x[i];
							float dy = y[i + 1] - y[i];
							float mag = sqrt(dx*dx + dy*dy);
							theta[i] = atan2f(dy / mag, dx / mag);
						}
						else
						{
							theta[i] = theta[i - 1];
						}

						//.. put into linked list
						int icell = int(floor(x[i] / dcell));
						int jcell = int(floor(y[i] / dcell));
						int scell = jcell*nxcell + icell;
						ptz[i] = heads[scell];
						heads[scell] = i;
					}
				}

				// Dump corner particles at end of file to trash *****************************

				for (int t = 0; t < 4; t++)
					fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;

				// Components of calculation *************************************************

				std::vector<float> G;
				G.resize(numBin + 1, 0.0f);
				float aveOver = 1;

				//.. calculate using linked list *********************************************

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

									//.. MAIN CALCULATION ************************************

									if (r2 <= intRange2)
									{
										//.. find bin
										float r = sqrt(r2);
										int b = int(r / binWidth);
										if ((b >= 0) && (b < numBin))
										{
											//.. relative angle between worms
											float dtheta = theta[ii] - theta[jj];

											//.. 2pi boundary conditons
											if (dtheta > _PI) dtheta -= _2PI;
											if (dtheta < -_PI) dtheta += _2PI;

											//.. calculate <n(r)*n(r+dr)>
											float vv = abs(cosf(dtheta));
											aveOver += 1;
											G[b] += vv;
										}
									}
									jj = ptz[jj];
								}
								ii = ptz[ii];
							}
						}
					}
				}

				//.. Average G(r) and print to out file *************************************

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

				//.. eliminate dynamic memory ***********************************************
				delete[] x;
				delete[] y;
				delete[] theta;
				delete[] ptz;
				delete[] heads;
			} // End of File

			fin.close();
		}
		fout.close();
		return EXIT_SUCCESS;
	}
}