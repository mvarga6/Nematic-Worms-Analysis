// Nematic Worm Analysis
// 11.11.14

/*
'follow' stand-alone analysis program
Track the positions of a list of particles into individual files.
*/

#define _USE_MATH_DEFINES

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
		<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
		<< "\t-e,--end\t\tDefine id of last input file (default=1)\n"
		<< "If no in/out options specified, default output file name is 'unnamed.txt'\n"
		<< "and input file must follow program instance\n"
		<< std::endl;
}

static int process_arg(std::string		&basename,
					   std::string		&foutbasename,
					   const int		argc,
					   char				*argv[],
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
	int		startFileId = 1;
	int		endFileId = 1;
	int		numFrame = 0;

	//.. process command line arguments
	process_arg(finBaseName, foutBaseName, argc, argv, startFileId, endFileId);

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

		//.. make output file for input file to be converted into
		std::ostringstream ofname;
		ofname << foutBaseName << fid << ".xyz";
		std::ofstream fout(ofname.str(), std::ios::out);

		//.. stuff
		int		numParticles;
		int		numPerWorm;
		int		numWorms;
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

			//.. data containers
			char  * c = new char[numParticles];
			float * x = new float[numParticles];
			float * y = new float[numParticles];

			//.. corner of box
			char  * bc = new char[4];
			float * bx = new float[4];
			float * by = new float[4];

			//.. center of frame position
			const float cofX = hx / 2.0f;
			const float cofY = hy / 2.0f;

			// Read in X,Y,Z positions for frame
			std::cout << "\nReading frame " << numFrame << " from " << ifname.str() << std::endl;
			int t = 0;
			for (int i = 0; i < numParticles; i++)
			{
				fin >> c[i] >> x[i] >> y[i] >> floatTrash;
  			}

			// Dump corner particles at end of file to trash
			for (int i = 0; i < 4; i++)
				fin >> bc[i] >> bx[i] >> by[i] >> floatTrash;

			// Calculate average position in mapped coords
			float gx_ave = 0.0f;
			float fx_ave = 0.0f;
			float gy_ave = 0.0f;
			float fy_ave = 0.0f;
			for (int i = 0; i < numParticles; i++)
			{
				float thetax = (x[i] / hx) * 2.0f * M_PI;
				gx_ave += cos(thetax);
				fx_ave += sin(thetax);

				float thetay = (y[i] / hy) * 2.0f * M_PI;
				gy_ave += cos(thetay);
				fy_ave += sin(thetay);
			}
			gx_ave /= float(numParticles);
			fx_ave /= float(numParticles);
			gy_ave /= float(numParticles);
			fy_ave /= float(numParticles);

			const float thetax_ave = atan2f(-fx_ave, -gx_ave) + M_PI;
			const float thetay_ave = atan2f(-fy_ave, -gy_ave) + M_PI;

			//.. convert back to cartesian coords
			float comX = (hx * thetax_ave) / (2.0f * M_PI);
			float comY = (hy * thetay_ave) / (2.0f * M_PI);

			//.. check that center of mass is in box
			if (comX < 0 || comX > hx || comY < 0 || comY > hy)
			{
				std::cout << "WARNING:\tCenter of mass outside of system bounds.\n";
				std::cout << "\t\tMoving through periodic boundary.\n";

				std::cout << "\t\tcomX = " << comX << " --> ";
				if (comX > hx)	 comX -= hx;
				if (comX < 0.0f) comX += hx;
				std::cout << comX << std::endl;

				std::cout << "\t\tcomY = " << comY << " --> ";
				if (comY > hy)	 comY -= hy;
				if (comY < 0.0f) comY += hy;
				std::cout << comY << std::endl;
			}

			//.. determine the necessary shift of center of frame
			const float shiftX = cofX - comX;
			const float shiftY = cofY - comY;

			//.. shift and apply periodic boundaries
			for (int i = 0; i < numParticles; i++)
			{
				x[i] += shiftX;
				y[i] += shiftY;

				if (x[i] > hx) x[i] -= hx;
				if (x[i] < 0.0f) x[i] += hx;
				if (y[i] > hy) y[i] -= hy;
				if (y[i] < 0.0f) y[i] += hy;
			}

			//.. print to output file
			fout << numParticles + 4 << std::endl;
			fout << numPerWorm << " " << k2spring << " " << hx << " " << hy << std::endl;
			for (int i = 0; i < numParticles; i++)
			{
				fout << c[i] << " " << x[i] << " " << y[i] << " 0\n";
			}

			//.. print box corners
			for (int i = 0; i < 4; i++)
			{
				fout << bc[i] << " " << bx[i] << " " << by[i] << " 0\n";
			}

			delete[] c;
			delete[] x;
			delete[] y;
			delete[] bc;
			delete[] bx;
			delete[] by;
			numFrame++;
		}
		fin.close();
		fout.close();
	}

	return EXIT_SUCCESS;
}