// Nematic Worm Analysis
// 11.11.14

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <math.h>

namespace svden {

	const float _PI = 3.14159265359f;
	const float _2PI = 2 * _PI;

	const std::string funcName = "svden";

	static void show_usage(std::string name)
	{
		std::cerr << "\a";
		std::cerr << "Usage: " << name << " <option(s)> file...\n"
			<< "\t-i,--input\t\tInput file name\n"
			<< "\t-o,--output\t\tOutput file name\n"
			<< "Options:\n"
			<< "\t-h,--help\t\tShow this help message\n"
			<< "\t-r,--range\t\tCalculation interaction range (default=5.0)\n"
			<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
			<< "\t-e,--end\t\tDefine id of last input file (default=1)\n"
			<< "If no in/out options specified, default output file name is 'Q-tensor.txt'\n"
			<< "and input file must follow program instance\n"
			<< std::endl;
	}

	static int process_arg(std::string		&basename,
		std::ofstream 	&fout,
		std::vector<std::string> argv,
		double 			&range,
		int				&sfid,
		int				&efid)
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
					basename = argv[i + 1];
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 2;
				}
			}
			else if ((arg == "-r") || (arg == "--range")){
				if (i + 1 < argc){
					range = strtod(argv[i + 1].c_str(), NULL);
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 3;
				}
			}
			else if ((arg == "-s") || (arg == "--start")){
				if (i + 1 < argc){
					sfid = int(strtod(argv[i + 1].c_str(), NULL));
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 4;
				}
			}
			else if ((arg == "-e") || (arg == "--end")){
				if (i + 1 < argc){
					efid = int(strtod(argv[i + 1].c_str(), NULL));
					i++;
				}
				else {
					std::cerr << "--input option requires one argument." << std::endl;
					return 5;
				}
			}
			else if ((arg == "-o") || (arg == "--output")){
				if (i + 1 < argc){
					fout.open(argv[i + 1]);
					i++;
				}
				else{
					std::cerr << "--output option requires one argument." << std::endl;
					return 6;
				}
			}
			else{
				fout.open(funcName + ".csv");
				return 0;
			}
		}
		return 0;
	}

	// Returns r2Com value
	/* Notes:   dx,dy,_dx,_dy must be initiallized with x1-x2.

	float fullBoundaryConditionCalculation(float&  x1, float&  x2, float&  y1, float&  y2, float&  dx, float&  dy,
	float& _x1, float& _x2, float& _y1, float& _y2, float& _dx, float& _dy,
	float& Xcom1, float& Ycom1,
	float& Xcom2, float& Ycom2
	)
	{
	if(dx > hxo2){
	dx -= hx;
	x1 -= hx;
	}
	if(dx <-hxo2){
	dx += hx;
	x1 += hx;
	}
	if(dy > hyo2){
	dy -= hy;
	y1 -= hy;
	}
	if(dy <-hyo2){
	dy += hy;
	y1 += hy;
	}
	// Periodic Boundaries
	float Xcom1 = (x1 + x2)/2.0;
	float Ycom1 = (y1 + y2)/2.0;
	if(Xcom1 > hx) Xcom1 -= hx;
	if(Xcom1 < 0)  Xcom1 += hx;
	if(Ycom1 > hy) Ycom1 -= hy;
	if(Ycom1 < 0)  Ycom1 += hy;
	}*/

	int calculate(std::vector<std::string> argv){

		//std::vector<std::ifstream*> fin;
		std::string finBaseName;
		std::ofstream fout;
		std::ofstream fout2("test_counts.csv");

		//.. user input parameters
		double	intRange = 5.0;
		int		startFileId = 1;
		int		endFileId = 1;
		int		numFrame = 1;
		float	binWidth = 0.001f;
		int		maxNumBin = 0;

		//.. process command line sub arguments
		int process_arg_status = process_arg(finBaseName, fout, argv, intRange, startFileId, endFileId);
		if (process_arg_status != 0) return process_arg_status;

		//.. linked list parameters (for worms COM, not particles)
		float dcell = float(intRange);
		int nxcell, nycell, ncell;
		const int ddx[5] = { 0, -1, 0, 1, 1 };
		const int ddy[5] = { 0, 1, 1, 1, 0 };
		int * ptz = { 0 };
		int * heads = { 0 };

		for (int fid = startFileId; fid <= endFileId; fid++)
		{
			//.. make ifstream and open
			std::ostringstream ifname;
			ifname << finBaseName << fid << ".xyz";
			std::ifstream fin(ifname.str(), std::ios::in);

			if (!fin.is_open()){
				std::cerr << "Error opening file '" << ifname.str() << "'\n"
					<< "Check for correct input name and/or directory\n" << std::endl;
				show_usage(funcName);
				return 20;
			}
			// Simulation parameters
			int		numParticles;
			int		numPerWorm;
			int		numWorms;
			double	numWorms2;
			char	charTrash;
			float	k2spring;
			float	hx;
			float	hy;
			float	floatTrash;
			float	intRange2 = intRange * intRange;
			float	intArea = _PI * intRange2;

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
				float hxo2 = hx / 2.0;
				float hyo2 = hy / 2.0;

				//.. make linked list
				nxcell = int(ceil(hx / dcell));
				nycell = int(ceil(hy / dcell));
				ncell = nxcell*nycell;

				// Allocate Containers
				float* comx = new float[numWorms];
				float* comy = new float[numWorms];
				float* theta = new float[numWorms];
				ptz = new int[numWorms];
				heads = new int[ncell];

				//.. init heads
				for (int i = 0; i < ncell; i++)
					heads[i] = -1;

				// Read in X,Y,Z positions for frame
				//.. and set linked list
				std::cout << "Reading frame " << numFrame << " from " << ifname.str() << std::endl;
				for (int w = 0; w < numWorms; w++)
				{
					//.. read in particles for worm w
					float * x = new float[numPerWorm];
					float * y = new float[numPerWorm];
					for (int p = 0; p < numPerWorm; p++)
						fin >> charTrash >> x[p] >> y[p] >> floatTrash;

					//.. get COM for w
					float x1 = x[0];
					float y1 = y[0];
					float x2 = x[numPerWorm - 1];
					float y2 = y[numPerWorm - 1];

					delete[] x;
					delete[] y;

					float dx = x1 - x2;
					float dy = y1 - y2;
					if (dx > hxo2){
						dx -= hx;
						x1 -= hx;
					}
					if (dx < -hxo2){
						dx += hx;
						x1 += hx;
					}
					if (dy > hyo2){
						dy -= hy;
						y1 -= hy;
					}
					if (dy < -hyo2){
						dy += hy;
						y1 += hy;
					}
					//.. calculate theta for worm w
					theta[w] = std::atan2(dy, dx);

					// Periodic Boundaries on COM
					float Xcom = (x1 + x2) / 2.0;
					float Ycom = (y1 + y2) / 2.0;
					if (Xcom > hx) Xcom -= hx;
					if (Xcom < 0)  Xcom += hx;
					if (Ycom > hy) Ycom -= hy;
					if (Ycom < 0)  Ycom += hy;

					//.. add to linked list and store COM
					int icell = int(floor(Xcom / dcell));
					int jcell = int(floor(Ycom / dcell));
					int scell = jcell*nxcell + icell;
					ptz[w] = heads[scell];
					heads[scell] = w;
					comx[w] = Xcom;
					comy[w] = Ycom;
				}

				// Dump corner particles at end of file to trash
				for (int t = 0; t < 4; t++)
					fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;

				// Components of calculation
				float totalOrderParameter = 0;
				float totalCos2Sum = 0;
				float totalCosSinSum = 0;
				float totalAveOver = 0;

				//.. histogram for frame
				std::vector<float> svd;
				std::vector<int> numInBin;

				//.. calculate using linked list
				std::cout << "Calculating ... \n";
				for (int ic = 0; ic < nxcell; ic++)
				{
					for (int jc = 0; jc < nycell; jc++)
					{
						//.. scalar add of cell
						int scell = jc*nxcell + ic;
						//std::cout << scell << "\t";

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

							int iw = heads[scell];
							while (iw >= 0)
							{
								// Components of calculation
								float orderParameter = 0;
								float cos2Sum = 0;
								float cosSinSum = 0;
								float aveOver = 0;

								int jw = heads[scnab];
								while (jw >= 0)
								{
									//std::cout << iw << ", " << jw << std::endl;

									// Boundary conditions on Xcoms
									float dXcom = comx[iw] - comx[jw];
									float dYcom = comy[iw] - comy[jw];
									if (dXcom > hxo2) dXcom -= hx;
									if (dXcom < -hxo2) dXcom += hx;
									if (dYcom > hyo2) dYcom -= hy;
									if (dYcom < -hyo2) dYcom += hy;

									//.. only if close enough
									const float R2com = dXcom*dXcom + dYcom*dYcom;
									if (R2com <= intRange2)
									{
										//.. relative angle between worms
										float dtheta = theta[iw] - theta[jw];

										//.. 2pi boundary conditons
										if (dtheta > _PI) dtheta -= _2PI;
										if (dtheta < -_PI) dtheta += _2PI;

										cos2Sum += std::cos(dtheta)*std::cos(dtheta);
										cosSinSum += std::cos(dtheta)*std::sin(dtheta);
										totalCos2Sum += cos2Sum;
										totalCosSinSum += cosSinSum;
										totalAveOver += 1;
										aveOver += 1;
									}
									jw = ptz[jw];
								}

								//.. add to histogram of S vs Density
								// Form Averages and calculate OrderParameter for Frame
								const float aveCos2 = cos2Sum / aveOver;
								const float aveCosSin = cosSinSum / aveOver;
								orderParameter = 2 * std::sqrt((aveCos2 - 0.5f)*(aveCos2 - 0.5f) + aveCosSin*aveCosSin);

								//.. find bin and resize if necessary
								int b = int((aveOver / intArea) / binWidth);
								//std::cout << b << "," << svd.size() << "\t";

								while (svd.size() <= b)
								{
									//std::cout << "in\t";
									svd.push_back(0.0f);
									numInBin.push_back(0);
								}

								//float stop;
								//std::cout << b << "," << svd.size() << "\t";
								//std::cin >> stop;

								//.. add S to bin b
								svd[b] += orderParameter;
								numInBin[b]++;

								//std::cin >> stop;

								//.. move to next particle
								iw = ptz[iw];
							}
						}
					}
				}
				/*
				//.. loop through all worms
				for (int w1 = 0; w1 < numWorms; w1++)
				{
				//.. get head to tail vector
				float x1 = x[w1][0];
				float y1 = y[w1][0];
				float x2 = x[w1][numPerWorm - 1];
				float y2 = y[w1][numPerWorm - 1];
				float dx = x1 - x2;
				float dy = y1 - y2;
				if (dx > hxo2){
				dx -= hx;
				x1 -= hx;
				}
				if (dx < -hxo2){
				dx += hx;
				x1 += hx;
				}
				if (dy > hyo2){
				dy -= hy;
				y1 -= hy;
				}
				if (dy < -hyo2){
				dy += hy;
				y1 += hy;
				}
				// Periodic Boundaries
				float Xcom1 = (x1 + x2) / 2.0;
				float Ycom1 = (y1 + y2) / 2.0;
				if (Xcom1 > hx) Xcom1 -= hx;
				if (Xcom1 < 0)  Xcom1 += hx;
				if (Ycom1 > hy) Ycom1 -= hy;
				if (Ycom1 < 0)  Ycom1 += hy;

				//.. loop through all other unconsidered pairs
				for (int w2 = w1 + 1; w2 < numWorms; w2++)
				{
				float _x1 = x[w2][0];
				float _x2 = x[w2][numPerWorm - 1];
				float _y1 = y[w2][0];
				float _y2 = y[w2][numPerWorm - 1];
				float _dx = _x1 - _x2;
				float _dy = _y1 - _y2;
				if (_dx > hxo2){
				_dx -= hx;
				_x1 -= hx;
				}
				if (_dx < -hxo2){
				_dx += hx;
				_x1 += hx;
				}
				if (_dy > hyo2){
				_dy -= hy;
				_y1 -= hy;
				}
				if (_dy < -hyo2){
				_dy += hy;
				_y1 += hy;
				}
				// Periodic Boundaries
				float Xcom2 = (_x1 + _x2) / 2.0;
				float Ycom2 = (_y1 + _y2) / 2.0;
				if (Xcom2 > hx) Xcom2 -= hx;
				if (Xcom2 < 0)  Xcom2 += hx;
				if (Ycom2 > hy) Ycom2 -= hy;
				if (Ycom2 < 0)  Ycom2 += hy;

				// Boundary conditions on Xcoms
				float dXcom = Xcom1 - Xcom2;
				float dYcom = Ycom1 - Ycom2;
				if (dXcom > hxo2) dXcom -= hx;
				if (dXcom < -hxo2) dXcom += hx;
				if (dYcom > hyo2) dYcom -= hy;
				if (dYcom < -hyo2) dYcom += hy;

				//.. check for close enough
				const float R2com = dXcom*dXcom + dYcom*dYcom;
				if (R2com > intRange2) continue;

				//.. update properties
				const float theta1 = std::atan2(dy, dx);
				const float theta2 = std::atan2(_dy, _dx);
				const float dtheta = theta1 - theta2;
				cos2Sum += std::cos(dtheta)*std::cos(dtheta);
				cosSinSum += std::cos(dtheta)*std::sin(dtheta);
				aveOver += 1;

				}
				}*/

				//.. print to file
				const int numBin = svd.size();
				fout2 << numFrame;
				fout << numFrame++;
				for (int i = 0; i < numBin; i++)
				{
					float printVal;
					if (numInBin.at(i) > 0)
					{
						printVal = svd[i] / float(numInBin.at(i));

					}
					else
					{
						printVal = 0.0f;
					}
					fout << ", " << printVal;
					fout2 << ", " << numInBin.at(i);
				}
				fout << std::endl;
				fout2 << std::endl;

				//.. store for printing bin labels
				if (numBin > maxNumBin) maxNumBin = numBin;

				//.. print to console
				totalOrderParameter = 2 * std::sqrtf(std::powf((totalCos2Sum / totalAveOver - 0.5f), 2.0f) + std::powf(totalCosSinSum / totalAveOver, 2.0f));
				std::cout << "S = " << totalOrderParameter << std::endl;
				fout.flush();
				//fout << numFrame++ << ", " << totalOrderParameter << std::endl;

				// Handle Dynamic Arrays
				delete[] comx;
				delete[] comy;
				delete[] theta;
				delete[] ptz;
				delete[] heads;
			} // End of File

			fin.close();

			//totalOrderParameter /= float(numFrame);
			//std::cout << "S_ave = " << totalOrderParameter << std::endl;

		}

		//.. print bin labels
		fout << "bins";
		for (int i = 1; i <= maxNumBin; i++)
			fout << ", " << i*binWidth;
		fout << std::endl;

		fout.close();
		fout2.close();
		return EXIT_SUCCESS;
	}
}