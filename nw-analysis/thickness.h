// Nematic Worms Analysis
// 2.4.16
// Detect thickness of worm layer
// Mike Varga
// -------------------------------
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <math.h>
// -------------------------------
namespace thickness_3d {
	const std::string funcName = "layer-thickness";
	
	// --------------------------------------------
	// Prints to user the program usage.
	static void show_usage(std::string name)
	{
		std::cerr << "\a";
		std::cerr << "Usage: " << name << " <option(s)> file...\n"
			<< "\t-i,--input\t\tAdd input file\n"
			<< "\t-o,--output\t\tOutput file name\n"
			<< "Options:\n"
			<< "\t-h,--help\t\tShow this help message\n"
			<< "\t-w,--boxwidth\t\tSet width of integration cells\n"
			<< "\t-d,--dimension\t\tSet width of integration cells\n"
			<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
			<< "\t-e,--end\t\tDefine id of last input file (default=1)\n"
			<< "If no in/out options specified, default output file name is 'unnamed.txt'\n"
			<< "and input file must follow program instance\n"
			<< std::endl;
	}
	
	// ------------------------------------------------
	// Process cmdline args and assign as needed.
	static int process_arg(std::string		&basename,
		std::ofstream	&fout,
		std::vector<std::string> argv,
		float			&boxwidth,
		int				&dim,
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
					fout.open(argv[i + 1]);
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
			else if ((arg == "-d") || (arg == "--dimension"))
			{
				if (i + 1 < argc){
					int assign = (int)strtof(argv[i + 1].c_str(), NULL);
					dim = assign;
					i++;
				}
				else
				{
					std::cerr << "--output option requires one argument." << std::endl;
					return 5;
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
				fout.open(funcName + ".csv");
				return 0;
			}
		}
		return 0;
	}

	// -------------------------------------------------
	// Run calculation
	int calculate(std::vector<std::string> argv){

		std::string finBaseName;
		std::ofstream fout;
		float * x0 = { 0 };
		float * y0 = { 0 };
		float	boxWidth = 2.0f;
		int		startFileId = 1;
		int		endFileId = 1;
		int		numFrame = 0;
		int		dim = 2; // Z by default

		//.. process and assign cmdline args
		int process_arg_status = process_arg(finBaseName, 
			fout, 
			argv, 
			boxWidth, 
			dim, 
			startFileId, 
			endFileId);
		if (process_arg_status != 0) return process_arg_status;

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
				show_usage(funcName);
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

				float * x = new float[numParticles];
				float * y = new float[numParticles];

				// Read in X,Y,Z positions for frame
				std::cout << "Reading frame " << numFrame << " from " << ifname.str() << std::endl;
				for (int i = 0; i < numParticles; i++)
					fin >> charTrash >> x[i] >> y[i] >> floatTrash;

				// Dump corner particles at end of file to trash
				for (int t = 0; t < 4; t++)
					fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;

				//.. alloc only at first frame
				if (numFrame == 0)
				{
					x0 = new float[numParticles];
					y0 = new float[numParticles];
					for (int i = 0; i < numParticles; i++)
					{
						x0[i] = x[i];
						y0[i] = y[i];
					}
				}

				// Calculate Average Properties
				float Vx = 0.0f;
				float Vy = 0.0f;
				for (int p = 0; p < numParticles; p++)
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
					Vx += dx / Dt;
					Vy += dy / Dt;
				}

				//.. print px py to file
				fout << Vx / numParticles << ", " << Vy / numParticles << std::endl;

				//.. store for next frame
				for (int i = 0; i < numParticles; i++)
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

		delete[] x0;
		delete[] y0;
		fout.close();

		return EXIT_SUCCESS;
	}
}