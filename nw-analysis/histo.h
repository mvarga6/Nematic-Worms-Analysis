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

namespace histo {

	const float _PI = 3.14159265359f;
	const float _2PI = 2 * _PI;

	const std::string funcName = "histo";

	static void show_usage(std::string name)
	{
		std::cerr << "\a";
		std::cerr << "Usage: " << name << " <option(s)> file...\n"
			<< "\t-i,--input\t\tInput file name\n"
			<< "\t-o,--output\t\tOutput file name\n"
			<< "Options:\n"
			<< "\t-h,--help\t\tShow this help message\n"
			<< "\t-N,--bins\t\tSplit histogram into N bins (default=45)\n"
			<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
			<< "\t-e,--end\t\tDefine id of last input file (default=1)\n"
			<< "If no in/out options specified, default output file name is 'hist.csv'\n"
			<< "and input file must follow program instance\n"
			<< std::endl;
	}

	static int process_arg(std::string		&basename,
		std::ofstream 	&fout,
		std::vector<std::string> argv,
		int 			&bins,
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
			else if ((arg == "-N") || (arg == "--bins")){
				if (i + 1 < argc){
					bins = int(strtod(argv[i + 1].c_str(), NULL));
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

	int calculate(std::vector<std::string> argv){

		//std::vector<std::ifstream*> fin;
		std::string finBaseName;
		std::ofstream fout;

		//.. user input parameters
		int		startFileId = 1;
		int		endFileId = 1;
		int     numFrame = 1;
		int		numBins = 45;

		//.. process command line sub arguments
		int process_arg_status = process_arg(finBaseName, fout, argv, numBins, startFileId, endFileId);
		if (process_arg_status != 0) return process_arg_status;

		const float binWidth = _2PI / float(numBins);
		fout << "bins";
		for (int i = 0; i < numBins; i++)
		{
			fout << "," << -_PI + i * binWidth;
		}
		fout << std::endl;

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

				//.. make histogram
				std::vector<int> histogram;
				histogram.resize(numBins, 0);

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
					float theta = std::atan2(dy, dx);

					// Periodic Boundaries on COM
					float Xcom = (x1 + x2) / 2.0;
					float Ycom = (y1 + y2) / 2.0;
					if (Xcom > hx) Xcom -= hx;
					if (Xcom < 0)  Xcom += hx;
					if (Ycom > hy) Ycom -= hy;
					if (Ycom < 0)  Ycom += hy;

					//.. add to histogram
					int b = int((theta + _PI) / binWidth);

					if (b < 0 || b >= histogram.size()) continue;

					histogram.at(b)++;
				}

				// Dump corner particles at end of file to trash
				for (int t = 0; t < 4; t++)
					fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;

				//.. print to output file
				fout << numFrame++;
				for (std::vector<int>::iterator it = histogram.begin();
					it != histogram.end(); ++it)
				{
					fout << "," << *it;
				}
				fout << std::endl;

			} // End of File

			fin.close();
		}
		fout.close();
		return EXIT_SUCCESS;
	}
}