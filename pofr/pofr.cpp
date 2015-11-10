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

const float _PI = 3.14159265359f;
const float _2PI = 2 * _PI;

static void show_usage(std::string name)
{
	std::cerr << "\a";
	std::cerr << "Usage: " << name << " <option(s)> file...\n"
		<< "\t-i,--input\t\tInput file name\n"
		<< "\t-o,--output\t\tOutput file name\n"
		<< "Options:\n"
		<< "\t-h,--help\t\tShow this help message\n"
		<< "\t-bw,--width\t\tBin width. dr, for P(r) (default=1.0)\n"
		<< "\t-n,--step\t\tDistribution after 'n' steps (default=50)\n"
		<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
		<< "\t-e,--end\t\tDefine id of last input file (default=1)\n"
		<< "If no in/out options specified, default output file name is 'Q-tensor.txt'\n"
		<< "and input file must follow program instance\n"
		<< std::endl;
}

static int process_arg(std::string		&basename,
	std::ofstream 	&fout,
	const int 		argc,
	char* 			argv[],
	float 			&width,
	int				&step,
	int				&sfid,
	int				&efid)
{
	for (int i = 1; i < argc; ++i){
		std::string arg = argv[i];
		if ((arg == "-h") || (arg == "--help")){
			show_usage(argv[0]);
			return 0;
		}
		else if ((arg == "-i") || (arg == "--input")){
			if (i + 1 < argc){
				//fin.push_back(new std::ifstream(argv[i + 1], std::ios::in));
				basename = argv[i + 1];
				i++;
			}
			else {
				std::cerr << "--input option requires one argument." << std::endl;
				return 1;
			}
		}
		else if ((arg == "-bw") || (arg == "--width")){
			if (i + 1 < argc){
				width = strtof(argv[i + 1], NULL);
				i++;
			}
			else {
				std::cerr << "--input option requires one argument." << std::endl;
				return 1;
			}
		}
		else if ((arg == "-n") || (arg == "--step")){
			if (i + 1 < argc){
				step = int(strtof(argv[i + 1], NULL));
				i++;
			}
			else {
				std::cerr << "--input option requires one argument." << std::endl;
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
		else if ((arg == "-o") || (arg == "--output")){
			if (i + 1 < argc){
				fout.open(argv[i + 1]);
				i++;
			}
			else{
				std::cerr << "--output option requires one argument." << std::endl;
				return 1;
			}
		}
		else{
			fout.open("SvT.csv");
			return 2;
		}
	}
}

/*************************************************************************************************
**************************************   MAIN FUNCTION   *****************************************
**************************************************************************************************/

int main(int argc, char* argv[]){

	//.. i/o things ************************************************************************
	std::string finBaseName;
	std::ofstream fout;

	//.. user input parameters *************************************************************
	float	binWidth = 1.0f;
	int		dstep = 10;
	int	   startFileId = 1;
	int	   endFileId = 1;
	int     numFrame = 0;
	process_arg(finBaseName, fout, argc, argv, binWidth, dstep, startFileId, endFileId);

	//.. number of bins and first line of data file ****************************************
	//const int numBin = int(intRange);
	//fout << "bins";
	//for (int i = 1; i <= numBin; i++)
	//	fout << ", " << i;
	//fout << std::endl;

	//.. Long lived containers *************************************************************
	//std::vector< float* > vx;
	//std::vector< float* > vy;
	//std::vector< float* > x_all;
	//std::vector< float* > y_all;
	std::vector< float* > dr_all;
	float * savex = { 0 };
	float * savey = { 0 };
	bool firstFrame = true;

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
			continue;
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
			//std::vector<float> x;
			//std::vector<float> y;
			//std::vector<float> framevx;
			//std::vector<float> framevy;
			float * x = new float[numParticles];
			float * y = new float[numParticles];
			//float * framedx = new float[numParticles];
			//float * framedy = new float[numParticles];
			float * framedr = new float[numParticles];

			// Read in X,Y,Z positions for frame and set linked list *********************
			std::cout << "Reading frame " << numFrame << " from " << ifname.str() << std::endl;
			for (int i = 0; i < numParticles; i++)
				fin >> charTrash >> x[i] >> y[i] >> floatTrash;

			// Dump corner particles at end of file to trash *****************************
			for (int t = 0; t < 4; t++)
				fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;

			//.. Just save positions the first time **************************************
			if (!firstFrame)
			{
				//.. calculate velocities ************************************************
				for (int i = 0; i < numParticles; i++)
				{
					float dx = savex[i] - x[i];
					float dy = savey[i] - y[i];
					if (dx > hxo2) dx -= hx;
					if (dx < -hxo2) dx += hx;
					if (dy > hyo2) dy -= hy;
					if (dy < -hxo2) dy += hy;

					//.. stored drs
					framedr[i] = sqrt(dx*dx + dy*dy);
					//framedx[i] = dx;
					//framedy[i] = dy;

					//.. save positions
					savex[i] = x[i];
					savey[i] = y[i];
				}
			}
			else
			{
				//.. alloc savex and savey
				savex = new float[numParticles];
				savey = new float[numParticles];

				//.. first frame saving then continue to load next frame *****************
				firstFrame = false;
				for (int i = 0; i < numParticles; i++)
				{
					savex[i] = x[i];
					savey[i] = y[i];
					//x_all.push_back(x);
					//y_all.push_back(y);
				}
				numFrame++;
				continue; // to next frame
			}

			//.. add to long lived container then remove first if full *******************
			if (dr_all.size() == dstep)
			{
				//.. add and erase
				//x_all.push_back(x);
				//y_all.push_back(y);
				dr_all.push_back(framedr);
				//x_all.erase(x_all.begin());
				//y_all.erase(y_all.begin());
				dr_all.erase(dr_all.begin());
			}
			else
			{
				//.. wait until we have sufficient
				//x_all.push_back(x);
				//y_all.push_back(y);
				dr_all.push_back(framedr);
				numFrame++;
				continue;
			}

			// Components of calculation *************************************************
			std::vector<float> P;
			float aveOver = 0;

			//.. Calculate G(t) *********************************************************
			std::cout << "Calculating...\n";
			for (int i = 0; i < numParticles; i++)
			{
				//.. calculate distanse
				float Ri = 0.0f;
				for (int j = 0; j < dr_all.size(); j++)
				{
					Ri += dr_all[j][i];
				}

				//.. get bin number
				const int bin = int(Ri / binWidth);

				//.. stop if negative somehow
				if (bin < 0) continue;

				//.. change size of P(r) if needed
				if (bin >= P.size())
				{
					int dSize = bin - P.size() + 1;
					for (int j = 0; j < dSize; j++)
					{
						P.push_back(0.0f);
					}
				}

				//.. add count to bin
				P.at(bin) += 1.0f;

				//.. add to G(t)
				aveOver += 1;
			}

			//.. Average P(r) and print to out file *************************************
			fout << numFrame++;
			for (int i = 0; i < P.size(); i++)
			{
				P[i] /= aveOver;
				fout << ", " << P[i];
			}
			fout << std::endl;

			delete[] x, y, framedr;

		} // End of File

		fin.close();
	}

	delete[] savex, savey;

	fout.close();
	return EXIT_SUCCESS;
}