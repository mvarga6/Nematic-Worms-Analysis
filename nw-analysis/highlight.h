// Nematic Worm Analysis
// 11.11.14

/*
'follow' stand-alone analysis program
Track the positions of a list of particles into individual files.
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <math.h>

namespace highlight {

	const std::string funcName = "highlight";

	static void show_usage(std::string name)
	{
		std::cerr << "\a";
		std::cerr << "Usage: " << name << " <option(s)> file...\n"
			<< "\t-i,--input\t\tAdd input file\n"
			<< "\t-o,--output\t\tOutput files base name\n"
			<< "Options:\n"
			<< "\t-h,--help\t\tShow this help message\n"
			<< "\t-t,--target\tAdd a target worm to label\n"
			<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
			<< "\t-e,--end\t\tDefine id of last input file (default=1)\n"
			<< "\t-l,--only\t\tPrint to output only the targets (default=false)\n"
			<< "If no in/out options specified, default output file name is 'unnamed.txt'\n"
			<< "and input file must follow program instance\n"
			<< std::endl;
	}

	static int process_arg(std::string		&basename,
		std::string		&foutbasename,
		std::vector<std::string> argv,
		std::vector<int> &targets,
		int				&sfid,
		int				&efid,
		bool			&solo)
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
			else if ((arg == "-o") || (arg == "--output"))
			{
				if (i + 1 < argc){
					foutbasename = argv[i + 1];
					i++;
				}
				else
				{
					std::cerr << "--output option requires one argument." << std::endl;
					return 3;
				}
			}
			else if ((arg == "-t") || (arg == "--target"))
			{
				if (i + 1 < argc){
					int assign = int(strtof(argv[i + 1].c_str(), NULL));
					targets.push_back(assign);
					i++;
				}
				else
				{
					std::cerr << "--output option requires one argument." << std::endl;
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
			else if ((arg == "-l") || (arg == "--only")){
				solo = true;
			}
			else
			{
				foutbasename = funcName;
				return 0;
			}
		}
		return 0;
	}

	bool allTrue(std::vector<bool> &vec)
	{
		for (int i = 0; i < vec.size(); i++)
		{
			if (vec[i] == false)
				return false;
		}
		return true;
	}

	bool checkInBounds(int id, int lower, int upper)
	{
		if (id < lower) return false;
		else if (id > upper) return false;
		else return true;
	}

	int calculate(std::vector<std::string> argv)
	{
		//.. file i/o names
		std::string finBaseName;
		std::string foutBaseName;

		//.. special stuff needed for entirety of calculation
		int	startFileId = 1;
		int	endFileId = 1;
		int	numFrame = 0;
		bool soloTargets = false;
		std::vector<int> targetList;

		//.. process command line arguments
		int process_arg_status = process_arg(finBaseName, foutBaseName, argv, targetList, startFileId, endFileId, soloTargets);
		if (process_arg_status != 0) return process_arg_status;

		//.. loop through all files
		const unsigned int numTargets = targetList.size();
		for (int fid = startFileId; fid <= endFileId; fid++)
		{
			//.. make sure it's open first
			std::ostringstream ifname;
			ifname << finBaseName << fid << ".xyz";
			std::ifstream fin(ifname.str(), std::ios::in);
			if (!fin.is_open())
			{
				std::cerr << "Error opening file" << std::endl
					<< "Check for correct input name and directory\n" << std::endl;
				show_usage(funcName);
				return 20;
			}

			//.. make output file for input file to be converted into
			std::ostringstream ofname;
			ofname << foutBaseName << fid << ".xyz";
			std::ofstream fout(ofname.str(), std::ios::out);

			//.. stuff
			int		numParticles;
			int		numPerWorm;
			int		numWorms;
			char	charTrash;
			float	k2spring;
			float	hx;
			float	hy;
			float	floatTrash;
			std::vector<bool> safeTarget;

			//.. each frame loop
			while (!fin.eof())
			{
				numParticles = 4; // Deals with case of black line at end of file.
				fin >> numParticles;
				fin >> numPerWorm >> k2spring >> hx >> hy;
				numParticles -= 4;

				if (numParticles == 0) break;

				numWorms = numParticles / numPerWorm;

				for (int i = 0; i < numTargets; i++)
				{
					safeTarget.push_back(checkInBounds(targetList[i], 0, numWorms - 1));
				}

				//.. check if targets are in correct range
				while (!allTrue(safeTarget))
				{
					//.. look through entire list
					for (int i = 0; i < numTargets; i++)
					{
						if (!safeTarget[i])
						{
							int newTarget = 0;
							std::cout << "\nWARNING:\tTarget worm out of range!\n";
							std::cout << "\t\tEnter new target in [0," << numWorms - 1 << "]:\t" << targetList[i] << " --> ";
							std::cin >> newTarget;
							targetList[i] = newTarget;
							safeTarget[i] = checkInBounds(newTarget, 0, numWorms - 1);
						}
					}
				}

				//.. particle positions and types
				float * x = new float[numParticles];
				float * y = new float[numParticles];

				//.. box corners
				float * bx = new float[4];
				float * by = new float[4];

				//.. read in X,Y,Z positions for frame
				std::cout << "\nReading frame " << numFrame << " from " << ifname.str() << std::endl;
				for (int i = 0; i < numParticles; i++)
				{
					fin >> charTrash >> x[i] >> y[i] >> floatTrash;
				}

				//.. dump corner particles at end of file to trash
				for (int t = 0; t < 4; t++)
					fin >> charTrash >> bx[t] >> by[t] >> floatTrash;


				//.. decide how to write output files
				if (!soloTargets) // normal case
				{
					fout << numParticles + 4 << std::endl;
					fout << numPerWorm << " " << k2spring << " " << hx << " " << hy << std::endl;
					for (int w = 0; w < numWorms; w++)
					{
						bool inList = false;
						for (int j = 0; j < numTargets; j++)
						{
							if (targetList[j] == w) inList = true;
						}

						//.. print entire worm
						char type = 'B';
						if (inList)
						{
							type = 'A';
						}

						//.. print worm with color "type"
						for (int p = 0; p < numPerWorm; p++)
						{
							int i = w*numPerWorm + p;
							fout << type << " " << x[i] << " " << y[i] << " 0" << std::endl;
						}
					}
				}
				else // solo targets case
				{
					fout << numTargets*numPerWorm + 4 << std::endl;
					fout << numPerWorm << " " << k2spring << " " << hx << " " << hy << std::endl;
					for (int t = 0; t < numTargets; t++)
					{
						int w = targetList[t];
						for (int p = 0; p < numPerWorm; p++)
						{
							int i = w*numPerWorm + p;
							fout << "A " << x[i] << " " << y[i] << " 0" << std::endl;
						}
					}
				}

				//.. print box corners
				for (int i = 0; i < 4; i++)
				{
					fout << "C " << bx[i] << " " << by[i] << " 0" << std::endl;
				}

				//.. delete containters
				delete[] x;
				delete[] y;
				delete[] bx;
				delete[] by;
				numFrame++;
			}
			fin.close();
			fout.close();
		}
		return EXIT_SUCCESS;
	}
}