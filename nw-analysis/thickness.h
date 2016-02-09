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
			<< "\t-i,--input\t\tAdd input file (w/o ext)\n"
			<< "\t-o,--output\t\tOutput file name (w/o ext)\n"
			<< "Options:\n"
			<< "\t-h,--help\t\tShow this help message\n"
			<< "\t-w,--boxwidth\t\tSet width of integration cells\n"
			<< "\t-d,--dimension\t\tSet width of integration cells\n"
			<< "\t-l,--layerpos\t\tPosition of layer to calculate thickness about\n"
			<< "\t-s,--start\t\tDefine id of first input file (default=1)\n"
			<< "\t-e,--end\t\tDefine id of last input file (default=1)\n"
			<< "If no in/out options specified, default output file name is 'unnamed.txt'\n"
			<< "and input file must follow program instance\n"
			<< std::endl;
	}
	
	// ------------------------------------------------
	// Process cmdline args and assign as needed.
	static int process_arg(std::string		&basename,
		std::string		&outbasename,
		std::vector<std::string> argv,
		float			&boxwidth,
		float			&layerpos,
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
					outbasename  = argv[i + 1];
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
			else if ((arg == "-l") || (arg == "--layerpos"))
			{
				if (i + 1 < argc){
					float assign = strtof(argv[i + 1].c_str(), NULL);
					layerpos = assign;
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
				outbasename = funcName;
				return 0;
			}
		}
		return 0;
	}

	// -------------------------------------------------
	// calculate root mean square displacement from layer
	// in uv-plane (rms of Z(u,v) )
	void rms_of_layer(float *u, float *v, 
		float *Z, const float& layerZ, 
		const int nparts, const float& boxwidth, 
		const int udim, const int vdim,
		float **result){
		
		//.. tmp memory for counting in each square
		int **cnts = new int*[udim];
		for (int i = 0; i < udim; i++){
			cnts[i] = new int[vdim];
			for (int j = 0; j < vdim; j++){
				cnts[i][j] = 1;
			}
		}

		//.. loop and calcuate r2 sum
		int b1, b2;
		for (int i = 0; i < nparts; i++){
			b1 = (int)(u[i] / boxwidth);
			b2 = (int)(v[i] / boxwidth);

			if (b1 < 0 || b2 < 0) continue;
			if (b1 >= udim || b2 >= vdim) continue;

			//.. calculate rms from layer
			result[b1][b2] += (Z[i] - layerZ)*(Z[i] - layerZ);
			cnts[b1][b2]++;
		}

		//.. turn result into sqrt of average in box
		for (int i = 0; i < udim; i++){
			for (int j = 0; j < vdim; j++){
				float ans = result[i][j] / float(cnts[i][j]);
				ans = std::sqrt(ans);
				result[i][j] = ans;
			}
		}

		//.. free tmp memory
		for (int i = 0; i < udim; i++){
			delete[] cnts[i];
		}
		delete[] cnts;
	}

	// -------------------------------------------------
	// Run calculation
	int calculate(std::vector<std::string> argv){

		std::string finBaseName;
		std::string foutBaseName;
		float * x = { 0 };
		float * y = { 0 };
		float * z = { 0 };
		float	boxWidth = 2.0f;
		float	layerPos = 0.0f;
		int		startFileId = -1;
		int		endFileId = -1;
		int		numFrame = 0;
		int		dim = 2; // Z by default

		//.. process and assign cmdline args
		int process_arg_status = process_arg(finBaseName, 
			foutBaseName, 
			argv, 
			boxWidth, 
			layerPos,
			dim, 
			startFileId, 
			endFileId);
		if (process_arg_status != 0) 
			return process_arg_status;

		//.. open output files
		std::ofstream fxyz(foutBaseName + ".xyzc", std::ios::out);
		std::ofstream fcsv(foutBaseName + ".csv", std::ios::out);
		if (!fxyz.is_open() || !fcsv.is_open())
			return 10;

		//.. loop through all files (or just once when -1 & -1)
		for (int fid = startFileId; fid <= endFileId; fid++)
		{
			//.. make sure it's open first
			std::ostringstream ifname;
			ifname << finBaseName;
			if (startFileId >= 0) ifname << fid;
			ifname << ".xyz";
			std::ifstream fin(ifname.str(), std::ios::in);
			if (!fin.is_open()){
				std::cerr << std::endl << funcName << ": Error opening file\n"
					<< funcName << ": Check for correct input name and directory";
				show_usage(funcName);
				return 20;
			}

			//.. stuff
			int		numParticles;
			int		numPerWorm;
			char	charTrash;
			float	k2spring;
			float	hx;
			float	hy;
			float	floatTrash;

			//.. each frame loop
			while (!fin.eof())
			{
				//.. read frame header
				std::string line;
				if (!std::getline(fin, line)) break;
				numParticles = (int)std::strtod(line.c_str(), NULL);
				printf("\n%s: Parts: %i", funcName.c_str(), numParticles);
				if (!std::getline(fin, line)) break;
				printf("\n\n%s: Comment line: %s", funcName.c_str(), line.c_str());
				numParticles -= 4;
				
				// allocate memory
				x = new float[numParticles];
				y = new float[numParticles];
				z = new float[numParticles];

				// Read in X,Y,Z positions for frame
				printf("\n%s: Reading frame %i from %s", funcName.c_str(), numFrame, ifname.str().c_str());
				for (int i = 0; i < numParticles; i++){
					std::getline(fin, line);
					std::stringstream row(line);
					row >> charTrash >> x[i] >> y[i] >> z[i];
				}

				// Dump 3 corner particles at end of file to trash
				for (int t = 0; t < 3; t++){
					std::getline(fin, line);
					std::stringstream row(line);
					row >> charTrash >> floatTrash >> floatTrash >> floatTrash;
				}

				// Get box size from fourth
				std::getline(fin, line);
				std::stringstream row(line);
				row >> charTrash >> hx >> hy >> floatTrash;

				// Make mean square displayment in desired dim
				const int xdim = (int)ceil(hx / boxWidth);
				const int ydim = (int)ceil(hy / boxWidth);
				float ** rmsqr_ = new float *[xdim];
				for (int i = 0; i < xdim; i++){
					rmsqr_[i] = new float[ydim];
					for (int j = 0; j < ydim; j++){
						rmsqr_[i][j] = 0.0f; // init to 0
					}
				}
				
				// Calculate root-mean-square displacement from layer position
				switch (dim) {

				case 0: // Slice of y-z plane
					rms_of_layer(y, z, x, layerPos, numParticles, boxWidth, ydim, 1, rmsqr_);
				case 1: // Slice of x-z plane
					rms_of_layer(x, z, y, layerPos, numParticles, boxWidth, xdim, 1, rmsqr_);
				case 2: // Slice of x-y place
					rms_of_layer(x, y, z, layerPos, numParticles, boxWidth, xdim, ydim, rmsqr_);
				}
				
				// Print to files
				fxyz << xdim * ydim << std::endl;
				fxyz << "Frame " << numFrame++ << std::endl;
				for (int j = 0; j < ydim; j++){
					for (int i = 0; i < xdim; i++){
						fxyz << "A " << i << " " << j << " 0 " << rmsqr_[i][j] << std::endl;
						fcsv << rmsqr_[i][j] << ", ";
					}
					fcsv << std::endl;
				}

				for (int i = 0; i < xdim; i++)
					delete[] rmsqr_[i];
				delete[] rmsqr_;
				printf("\n%s: Frame memory deleted.", funcName.c_str());
			}
			fin.close();
			printf("\n%s: Input file closed.", funcName.c_str());
		}
		fxyz.close();
		fcsv.close();
		printf("\n%s: Output files closed.", funcName.c_str());
		return EXIT_SUCCESS;
	}
}