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

namespace gtcorr {

	const float _PI = 3.14159265359f;
	const float _2PI = 2 * _PI;

	const std::string funcName = "gtcorr";

	static void show_usage(std::string name)
	{
		std::cerr << "\a";
		std::cerr << "Usage: " << name << " <option(s)> file...\n"
			<< "\t-i,--input\t\tInput file name\n"
			<< "\t-o,--output\t\tOutput file name\n"
			<< "Options:\n"
			<< "\t-h,--help\t\tShow this help message\n"
			<< "\t-m,--max\t\tMaximum correlation time (default=10.0)\n"
			<< "\t-N,--targets\t\tNumber of random targets (default=1000)\n"
			<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
			<< "\t-e,--start\t\tDefine id of last input file (default=1)\n"
			<< "If no in/out options specified, default output file name is 'Q-tensor.txt'\n"
			<< "and input file must follow program instance\n"
			<< std::endl;
	}

	static int process_arg(std::string		&basename,
		std::ofstream 	&fout,
		std::vector<std::string> argv,
		float 			&range,
		int				&sfid,
		int				&efid,
		int				&num)
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
			else if ((arg == "-N") || (arg == "--targets")){
				if (i + 1 < argc){
					num = int(strtod(argv[i + 1].c_str(), NULL));
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
		int	startFileId = 1;
		int	endFileId = 1;
		int numFrame = 1;
		int	numTarget = 1000;

		//.. process command line sub arguments
		int process_arg_status = process_arg(finBaseName, fout, argv, intRange, startFileId, endFileId, numTarget);
		if (process_arg_status != 0) return process_arg_status;

		//.. number of bins and first line of data file ****************************************

		const int numBin = int(intRange);
		fout << "bins";
		for (int i = 1; i <= numBin; i++)
			fout << ", " << i;
		fout << std::endl;

		//.. loop through all input files ******************************************************

		std::vector< std::vector<float> > theta;
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

				// Allocate Containers *******************************************************

				float * x = new float[numParticles];
				float * y = new float[numParticles];
				std::vector<float> frameTheta;

				// Read in X,Y,Z positions for frame and set linked list *********************
				std::cout << "Reading frame " << numFrame << " from " << ifname.str() << std::endl;
				for (int w = 0; w < numWorms; w++)
				{
					//.. read in postions for worm w *****************************************

					for (int p = 0; p < numPerWorm; p++)
					{
						int i = w*numPerWorm + p;
						fin >> charTrash >> x[i] >> y[i] >> floatTrash;
					}

					//.. calculate theta *****************************************************

					for (int p = 0; p < numPerWorm; p++)
					{
						int i = w*numPerWorm + p;
						if (p < numParticles - 1)
						{
							float dx = x[i + 1] - x[i];
							float dy = y[i + 1] - y[i];
							float mag = sqrt(dx*dx + dy*dy);
							frameTheta.push_back(atan2f(dy / mag, dx / mag));
						}
						else
						{
							float back = frameTheta.back();
							frameTheta.push_back(back);
						}
					}
				}

				// Dump corner particles at end of file to trash *****************************

				for (int t = 0; t < 4; t++)
					fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;

				//.. figure what to do with frameTheta ***************************************

				if (theta.size() == numBin)
				{
					//.. add and erase
					theta.push_back(frameTheta);
					theta.erase(theta.begin());
				}
				else
				{
					//.. just add
					theta.push_back(frameTheta);
				}

				// Components of calculation *************************************************

				std::vector<float> G;
				G.resize(theta.size(), 0.0f);
				float aveOver = 0;

				//.. Calculate G(t) *********************************************************
				std::cout << "Calculating...\n";
				for (int dt = 0; dt < G.size(); dt++)
				{
					for (int i = 0; i < numParticles; i++)
					{
						//.. change in angle
						float dtheta = theta[dt][i] - theta[0][i];

						//.. 2pi boundary conditons
						if (dtheta > _PI) dtheta -= _2PI;
						if (dtheta < -_PI) dtheta += _2PI;

						//.. add to G(t)
						aveOver += 1;
						G[dt] += cosf(dtheta);
					}
				}

				//.. Average G(t) and print to out file *************************************

				fout << numFrame++;
				for (int i = 0; i < G.size(); i++)
				{
					G[i] /= aveOver;
					fout << ", " << G[i];
				}
				fout << std::endl;

				//.. eliminate dynamic memory ***********************************************
				delete[] x;
				delete[] y;
			} // End of File

			fin.close();
		}
		fout.close();
		return EXIT_SUCCESS;
	}
}