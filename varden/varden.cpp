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
		<< "\t-l,--lower\t\tLower box size limit(default=1.0)\n"
		<< "\t-u,--upper\t\tUpper box size limit(default=25.0)\n"
		<< "\t-N,--boxes\t\tNumber of box sizes(default=25)\n"
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
	float 			&min,
	float			&max,
	int				&num,
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
		else if ((arg == "-l") || (arg == "--lower")){
			if (i + 1 < argc){
				min = strtof(argv[i + 1], NULL);
				i++;
			}
			else {
				std::cerr << "--input option requires one argument." << std::endl;
				return 1;
			}
		}
		else if ((arg == "-u") || (arg == "--upper")){
			if (i + 1 < argc){
				max = strtof(argv[i + 1], NULL);
				i++;
			}
			else {
				std::cerr << "--input option requires one argument." << std::endl;
				return 1;
			}
		}
		else if ((arg == "-N") || (arg == "--boxes")){
			if (i + 1 < argc){
				num = int(strtof(argv[i + 1], NULL));
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

	int		startFileId = 1;
	int		endFileId	= 1;
	int     numFrame	= 0;
	int		numBoxSize	= 25;
	float	minBoxSize	= 1.0f;
	float	maxBoxSize	= 25.0f;
	process_arg(finBaseName, fout, argc, argv, minBoxSize, maxBoxSize, numBoxSize, startFileId, endFileId);

	//.. calc box sizes ******************************************************************

	float dSize = (maxBoxSize - minBoxSize) / float(numBoxSize);

	std::vector<float> boxSize;
	boxSize.push_back(minBoxSize);
	std::cout << "\tBox Sizes:\n";
	for (int i = 1; i < numBoxSize; i++)
	{
		float prevSize = boxSize.back();
		boxSize.push_back(prevSize + dSize);
		std::cout << "\t" << boxSize.back();
	}
	std::cout << std::endl << std::endl;

	//.. total average stuff *************************************************************

	std::vector<float> totalDN;
	std::vector<float> totalAveN;
	totalDN.resize(numBoxSize, 0.0f);
	totalAveN.resize(numBoxSize, 0.0f);

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

			// Read in X,Y,Z positions for frame and set linked list *********************
			std::cout << "Reading frame " << numFrame++ << " from " << ifname.str() << std::endl;
			for (int i = 0; i < numParticles; i++)
				fin >> charTrash >> x[i] >> y[i] >> floatTrash;
			
			// Dump corner particles at end of file to trash *****************************

			for (int t = 0; t < 4; t++)
				fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;

			// Components of calculation *************************************************

			std::vector<float> DN;
			std::vector<float> aveN;

			//.. calculate ***************************************************************

			for (int i = 0; i < numBoxSize; i++)
			{
				int xbox = int(hx / boxSize[i]);
				int ybox = int(hy / boxSize[i]);
				int nbox = xbox*ybox;
				int * counts = new int[nbox];

				for (int b = 0; b < nbox; b++)
					counts[b] = 0;

				for (int j = 0; j < numParticles; j++)
				{
					int ibox = int(x[j] / boxSize[i]);
					int jbox = int(y[j] / boxSize[i]);
					if (ibox >= xbox || jbox >= ybox) continue;

					int sbox = jbox*xbox + ibox;
					if (sbox >= 0 && sbox < nbox)
					{
						counts[sbox]++;
					}
				}

				//.. calculate sum
				int sum = 0;
				for (int b = 0; b < nbox; b++)
				{
					sum += counts[b];
				}

				float nbar = float(sum) / float(nbox);
				aveN.push_back(nbar);

				float sqrtSum = 0.0f;
				for (int b = 0; b < nbox; b++)
				{
					float dn = float(counts[b]) - nbar;
					sqrtSum += dn*dn;
				}

				DN.push_back(sqrtf(sqrtSum / float(nbox - 1)));
			}
			
			//.. make printable vector
			const int numOfXValues = 25;
			//std::vector<float> xValues;
			std::vector<float> yValues;
			std::vector<float> nums;
			
			//xValues.resize(numOfXValues, 0.0f);
			yValues.resize(numOfXValues, 0.0f);
			nums.resize(numOfXValues, 0.0f);

			//.. get max and min N
			float NMax = 0;
			float NMin = 100;
			for (int j = 0; j < aveN.size(); j++)
			{
				if (aveN[j] > NMax) NMax = aveN[j];
				if (aveN[j] < NMin) NMin = aveN[j];
			}

			float dx = (NMax - NMin) / float(numOfXValues);
			for (int j = 0; j < aveN.size(); j++)
			{
				int xbin = int(aveN[j] / dx);
				if (xbin >= numOfXValues || xbin < 0) continue;

				nums[xbin] += 1.0f;
				yValues[xbin] += DN[j];
			}

			fout << "0";
			for (int j = 1; j < numOfXValues; j++)
			{
				fout << "," << j*dx;
			}
			fout << std::endl;

			if (nums[0] != 0.0f) fout << yValues[0] / nums[0];
			else fout << 0.0f;
			for (int j = 1; j < numOfXValues; j++)
			{
				if (nums[j] != 0.0f) fout << "," << yValues[j] / nums[j];
				else fout << 0.0f;
			}
			fout << std::endl << std::endl;

			//.. eliminate dynamic memory ************************************************
			delete[] x;
			delete[] y;
		} // End of File

		fin.close();
	}
	fout.close();
	return EXIT_SUCCESS;
}