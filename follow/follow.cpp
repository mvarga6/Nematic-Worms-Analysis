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

static void show_usage(std::string name)
{
	std::cerr << "\a";
	std::cerr << "Usage: " << name << " <option(s)> file...\n"
		<< "\t-i,--input\t\tAdd input file\n"
		<< "\t-o,--output\t\tOutput files base name\n"
		<< "Options:\n"
		<< "\t-h,--help\t\tShow this help message\n"
		<< "\t-t,--target\tAdd a target particle to follow\n"
		<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
		<< "\t-e,--end\t\tDefine id of last input file (default=1)\n"
		<< "If no in/out options specified, default output file name is 'unnamed.txt'\n"
		<< "and input file must follow program instance\n"
		<< std::endl;
}

static int process_arg(std::string		&basename,
	std::string		&foutbasename,
	const int		argc,
	char			*argv[],
	std::vector<int> &targets,
	int				&sfid,
	int				&efid)
{
	for (int i = 1; i < argc; ++i)
	{
		std::string arg = argv[i];
		if ((arg == "-h") || (arg == "--help"))
		{
			show_usage(argv[0]);
			return 0;
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
				return 1;
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
				return 1;
			}
		}
		else if ((arg == "-t") || (arg == "--target"))
		{
			if (i + 1 < argc){
				int assign = int(strtof(argv[i + 1], NULL));
				targets.push_back(assign);
				i++;
			}
			else
			{
				std::cerr << "--output option requires one argument." << std::endl;
				return 1;
			}
		}
		else if ((arg == "-s") || (arg == "--start")){
			if (i + 1 < argc){
				sfid = int(strtod(argv[i + 1], NULL));
				i++;
			}
			else {
				std::cerr << "--input option requires one argument." << std::endl;
				return 1;
			}
		}
		else if ((arg == "-e") || (arg == "--end")){
			if (i + 1 < argc){
				efid = int(strtod(argv[i + 1], NULL));
				i++;
			}
			else {
				std::cerr << "--input option requires one argument." << std::endl;
				return 1;
			}
		}
		else
		{
			foutbasename = "follow";
			return 2;
		}
	}
}

int main(int argc, char* argv[])
{
	//.. file i/o names
	std::string finBaseName;
	std::string foutBaseName;

	//.. special stuff needed for entirety of calculation
	float * x0 = { 0 };
	float * y0 = { 0 };
	float * xTot = { 0 };
	float * yTot = { 0 };
	int		startFileId = 1;
	int		endFileId = 1;
	int		numFrame = 0;
	std::vector<int> targetList;

	//.. process command line arguments
	process_arg(finBaseName, foutBaseName, argc, argv, targetList, startFileId, endFileId);

	//.. create output file for each target
	std::vector<std::ofstream*> fout;
	const unsigned int numTargets = targetList.size();
	for (int i = 0; i < numTargets; i++)
	{
		std::ostringstream ofname;
		ofname << foutBaseName << targetList[i] << ".csv";
		fout.push_back(new std::ofstream(ofname.str(), std::ios::out));
	}

	//.. confirm sizes match
	if (numTargets != fout.size())
	{
		std::cerr << "Error linking targets to files\n";
		return 10;
	}

	//.. loop through all files
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
			show_usage(argv[0]);
			return 20;
		}

		//.. stuff
		int		numParticles;
		int		numPerWorm;
		int		numWorms;
		double	numWorms2;
		char	charTrash;
		float	k2spring;
		float	hx;
		float	hy;
		float	floatTrash;

		//.. each frame loop
		while (!fin.eof())
		{
			numParticles = 4; // Deals with case of black line at end of file.
			fin >> numParticles;
			fin >> numPerWorm >> k2spring >> hx >> hy;
			numParticles -= 4;

			if (numParticles == 0) break;

			numWorms = numParticles / numPerWorm;
			numWorms2 = numWorms*numWorms;

			float * x = new float[numTargets];
			float * y = new float[numTargets];

			// Read in X,Y,Z positions for frame
			std::cout << "\nReading frame " << numFrame << " from " << ifname.str() << std::endl;
			int t = 0;
			for (int i = 0; i < numParticles; i++)
			{
				//.. is i in targets list?
				bool found = false;
				for (int c = 0; c < numTargets; c++)
				{
					if (targetList[c] == i) found = true;
				}

				if (found)
				{
					fin >> charTrash >> x[t] >> y[t] >> floatTrash;
					std::cout << "\nTarget # " << t << " found at { " << x[t] << "," << y[t] << " }";
					t++;
				}
				else
					fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;
			}

			// Dump corner particles at end of file to trash
			for (int t = 0; t < 4; t++)
				fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;

			//.. alloc only at first frame
			if (numFrame == 0)
			{
				x0 = new float[numTargets];
				y0 = new float[numTargets];
				xTot = new float[numTargets];
				yTot = new float[numTargets];
				for (int i = 0; i < numTargets; i++)
				{
					x0[i] = x[i];
					y0[i] = y[i];
					xTot[i] = 0.0f;
					yTot[i] = 0.0f;
  				}
			}

			// Calculate Average Properties
			for (int p = 0; p < numTargets; p++)
			{
				//.. raw distance travelled
				float dx = x[p] - x0[p];
				float dy = y[p] - y0[p];

				//.. boundary conditions
				if (dx > hx / 2.0f) dx -= hx;
				if (dx < -hx / 2.0f) dx += hx;
				if (dy > hy / 2.0f) dy -= hy;
				if (dy < -hy / 2.0f) dy += hy;

				//.. add to cummulative distance vector
				xTot[p] += dx;
				yTot[p] += dy;
			}

			//.. print to output files
			for (int i = 0; i < numTargets; i++)
			{
				*(fout[i]) << xTot[i] << ", " << yTot[i] << std::endl;
			}

			//.. store for next frame
			for (int i = 0; i < numTargets; i++)
			{
				x0[i] = x[i];
				y0[i] = y[i];
			}

			delete[] x;
			delete[] y;
			numFrame++;
		}
		fin.close();
	}

	//.. delete containers
	delete[] x0;
	delete[] y0;
	delete[] xTot;
	delete[] yTot;
	for (int i = 0; i < numTargets; i++)
	{
		fout[i]->close();
		delete fout[i];
	}

	return EXIT_SUCCESS;
}