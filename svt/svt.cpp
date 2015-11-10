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
		<< "\t-r,--range\t\tCalculation interaction range (default=5.0)\n"
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
	double 			&range,
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
		else if ((arg == "-r") || (arg == "--range")){
			if (i + 1 < argc){
				range = strtod(argv[i + 1], NULL);
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

int main(int argc, char* argv[]){

	//std::vector<std::ifstream*> fin;
	std::string finBaseName;
	std::ofstream fout;

	//.. user input parameters
	double intRange = 5.0;
	int	   startFileId = 1;
	int	   endFileId = 1;
	int     numFrame = 1;
	float	totalOrderParameter = 0;
	process_arg(finBaseName, fout, argc, argv, intRange, startFileId, endFileId);

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
			continue;
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
			float orderParameter = 0;
			float cos2Sum = 0;
			float cosSinSum = 0;
			float aveOver = 0;

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
									aveOver += 1;
								}

								jw = ptz[jw];
							}
							iw = ptz[iw];
						}
					}
				}
			}
			
			// Form Averages and calculate OrderParameter for Frame
			const float aveCos2 = cos2Sum / aveOver;
			const float aveCosSin = cosSinSum / aveOver;
			orderParameter = 2 * std::sqrt((aveCos2 - 0.5)*(aveCos2 - 0.5) + aveCosSin*aveCosSin);
			totalOrderParameter += orderParameter;

			std::cout << "S = " << orderParameter << std::endl;
			fout.flush();
			fout << numFrame++ << ", " << orderParameter << std::endl;

			// Handle Dynamic Arrays
			delete[] comx;
			delete[] comy;
			delete[] theta;
			delete[] ptz;
			delete[] heads;
		} // End of File

		fin.close();

		totalOrderParameter /= float(numFrame);
		std::cout << "S_ave = " << totalOrderParameter << std::endl;

	}
	fout.close();
	return EXIT_SUCCESS;
}