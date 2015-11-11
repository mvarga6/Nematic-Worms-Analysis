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
		<< "\t-bw,--width\t\tBin width for histogram (default=0.25)\n"
		<< "\t-d,--spacer\t\tParticle index spacing for Kuhn calc (default=2)\n"
		<< "\t-T,--kBT\t\tTemperature of simulation (default=1.0)\n"
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
	float			&width,
	float			&kBT,
	int				&spacer,
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
				width = float(strtod(argv[i + 1], NULL));
				i++;
			}
			else {
				std::cerr << "--input option requires one argument." << std::endl;
				return 1;
			}
		}
		else if ((arg == "-d") || (arg == "--spacer")){
			if (i + 1 < argc){
				spacer = int(strtod(argv[i + 1], NULL));
				i++;
			}
			else {
				std::cerr << "--input option requires one argument." << std::endl;
				return 1;
			}
		}
		else if ((arg == "-T") || (arg == "--kBT")){
			if (i + 1 < argc){
				kBT = strtof(argv[i + 1], NULL);
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
			fout.open("persist.csv");
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

	int	startFileId = 1;
	int	endFileId = 1;
	int numFrame = 1;
	int kuhnSpacer = 2;
	float binWidth = 0.25f;
	float kBT = 1.0f;
	process_arg(finBaseName, fout, argc, argv, binWidth, kBT, kuhnSpacer, startFileId, endFileId);

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

			float * x = new float[numParticles];
			float * y = new float[numParticles];
			float * bx = new float[numParticles];
			float * by = new float[numParticles];
			float * lp = new float[numWorms];

			// Read in X,Y,Z positions for frame and set linked list *********************
			std::cout << "Reading frame " << numFrame << " from " << ifname.str() << std::endl;
			//int i, j; // used for all particle indices (vs. worm and interworm indices)
			for (int w = 0; w < numWorms; w++)
			{
				//.. read in postions for worm w *****************************************

				for (int p = 0; p < numPerWorm; p++)
				{
					int i = w*numPerWorm + p;
					fin >> charTrash >> x[i] >> y[i] >> floatTrash;
				}

				//.. calculate bx and by *************************************************
				//float dx, dy;
				for (int p = 0; p < numPerWorm; p++)
				{
					int j = w*numPerWorm + p;
					if (p < numPerWorm - 1)
					{
						//.. Displacement
						float dx = x[j + 1] - x[j];
						float dy = y[j + 1] - y[j];
						
						//.. PBC
						if (dx > hxo2) dx -= hx;
						if (dx < -hxo2) dx += hx;
						if (dy > hyo2) dy -= hy;
						if (dy < -hyo2) dy += hy;

						//.. store bx and by
						bx[j] = dx; by[j] = dy;
					}
					else
					{
						bx[j] = bx[j - 1];
						by[j] = by[j - 1];
					}

					//std::cout << p << ": " << bx[j] << ", " << by[j] << std::endl;
				}
			}

			// Dump corner particles at end of file to trash *****************************

			for (int t = 0; t < 4; t++)
				fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;

			//.. calculate persistence lengths *********************************************

			std::cout << "Calculating ... \n";
			float ave_lp = 0.0f, max_lp = 0.0f, min_lp = 0.0f;
			for (int w = 0; w < numWorms; w++)
			{
				lp[w] = 0.0f;
				for (int p = 0; p < numPerWorm - 2; p++)
				{
					int i = numPerWorm*w + p;

					//.. dot product
					lp[w] += bx[i] * bx[i + 1] + by[i] * by[i + 1];
				}

				//.. determine max and min
				if (lp[w] < min_lp) min_lp = lp[w];
				if (lp[w] > max_lp) max_lp = lp[w];

				//.. add to average sum
				ave_lp += lp[w];
			}
			ave_lp /= float(numWorms);

			//.. calculate kuhn length
			/*float dx = 0.0f, dy = 0.0f;
			float tot_ave_P = 0.0f, dbl_dist_dot, ave_dist, ave_distOver;
			for (int w = 0; w < numWorms; w++)
			{
				//std::cout << std::endl << "WORM " << w << std::endl;
				//.. for kuhn length calculation
				ave_dist = 0.0f; ave_distOver = 0.0f; dbl_dist_dot = 0.0f;
				for (int p = 0; p < numPerWorm - kuhnSpacer; p++)
				{
					int i = numPerWorm*w + p;
					int ips = i + kuhnSpacer;

					dx = x[ips] - x[i];
					dy = y[ips] - y[i];

					//.. PBC
					if (dx > hxo2) dx -= hx;
					if (dx < -hxo2) dx += hx;
					if (dy > hyo2) dy -= hy;
					if (dy < -hyo2) dy += hy;

					ave_dist += sqrtf(dx*dx + dy*dy);
					dbl_dist_dot += bx[i] * bx[ips] + by[i] * by[ips];
					ave_distOver += 1.0f;
				}

				//.. calculate persistence length
				tot_ave_P += -(ave_dist / ave_distOver) / logf(dbl_dist_dot / ave_distOver);
			}
			tot_ave_P /= float(numWorms);*/

			std::cout << "\tKuhn length = " << ave_lp * 2.0f << std::endl
				<< "\tBending Stiffness = " << ave_lp * kBT << std::endl;
			
			//.. Histogram distribution from ave_lp *************************************

			std::cout << "\tProducing histogram\n";
			std::vector<float> hist;
			if (min_lp < 0) min_lp = 0.0f;
			int bin;
			float aveOver = 0.0f;
			for (int w = 0; w < numWorms; w++)
			{
				//.. find bin
				bin = int(lp[w] / binWidth);

				//.. fix possible negatives
				if (bin < 0) bin = 0;

				//.. adjust size if necessary
				if (bin >= hist.size())
				{
					int dsize = bin - hist.size() + 1;
					for (int add = 0; add < dsize; add++)
						hist.push_back(0.0f);
				}

				//.. add to bin
				hist.at(bin) += 1.0f;
				aveOver += 1.0f;
			}
			
			//.. print to file
			std::cout << "\tPrinting to data file\n";
			fout << numFrame++ << ", " << ave_lp << ", " << ave_lp * kBT;
			for (int c = 0; c < hist.size(); c++)
				fout << ", " << hist.at(c) / aveOver;
			fout << std::endl;

			//.. eliminate dynamic memory ***********************************************
			delete[] x, y;
			delete[] bx, by, lp;
		} // End of File

		fin.close();
	}
	fout.close();
	return EXIT_SUCCESS;
}
