// Nematic Worm Analysis
// 11.11.14

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <stdio.h>
#include <math.h>

class n {
public:

	n(void){
		nx  = 0;
		ny  = 0;
		mag = 0;
	}
	n(float x, float y){
		nx = x;
		ny = y;
		mag = std::sqrt(nx*nx + ny*ny);
		cosT = std::cos(atan(ny/nx));
	}
				
	void set(float x, float y){
		nx = x;
		ny = y;
		mag = std::sqrt(nx*nx + ny*ny);
		cosT = std::cos(atan(ny/nx));
	}
	
	float& getnx(void){ return nx;}
	float& getny(void){ return ny;}
	
	float operator*(const n& rhs){
		return (this->nx*rhs.nx + this->ny*rhs.ny);
	}
	
private:
	float nx;
	float ny;
	float mag;
	float cosT;
};

static void show_usage(std::string name)
{
	std::cerr << "\a";
    std::cerr << "Usage: " << name << " <option(s)> file...\n"
			  << "\t-i,--input\t\tInput file name\n"
			  << "\t-o,--output\t\tOutput file name\n"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
			  << "If no in/out options specified, default output file name is 'Q-tensor.txt'\n"
			  << "and input file must follow program instance\n"
              << std::endl;
}

static int process_arg(std::ifstream &fin, std::ofstream &fout, 
						const int argc, char* argv[])
{
	for (int i = 1; i < argc; ++i){
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")){
            show_usage(argv[0]);
            return 0;
        } 
		else if ((arg == "-i") || (arg == "--input")){
            if (i + 1 < argc){
                fin.open(argv[i + 1]);
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
            fin.open(argv[i]);
			fout.open("S.xyzc");
			return 2;
        }
    }
}

int main(int argc, char* argv[]){
	
	if(argc > 7){
		show_usage(argv[0]);
		return 1;
	}

	std::ifstream fin;
	std::ofstream fout;
	
	process_arg(fin,fout,argc,argv);

	if(!fin.is_open()){
		std::cerr << "Error opening file" << std::endl
				  << "Check for correct input name and directory\n" << std::endl;
		show_usage(argv[0]);
	}
	else{
		int		numParticles;
		int		numPerWorm; 
		int		numWorms;
		double	numWorms2;
		int     numFrame = 1;
		char	charTrash;
		float	k2spring;
		float	hx;
		float	hy;
		float	floatTrash;
		float	boxSize = 2.5f;
		
		while(!fin.eof()){
			numParticles = 4; // Deals with case of black line at end of file.
			fin >> numParticles;
			fin >> numPerWorm >> k2spring >> hx >> hy;
			numParticles -= 4;
			
			if(numParticles!=0){
			
				numWorms = numParticles/numPerWorm;
				numWorms2 = numWorms*numWorms;
			
				// Set up dynamic arrays 4 particles
				float** x = new float*[numWorms];
				float** y = new float*[numWorms];
				for(int w = 0; w < numWorms; w++){
					x[w] = new float[numPerWorm];
					y[w] = new float[numPerWorm];
				}
			
				// Determine lattice for calculating average S
				const int numXBox = int(std::ceil(hx/boxSize));
				const int numYBox = int(std::ceil(hy/boxSize));
				const int numBox  = numXBox * numYBox;
			
				// set up dynamic array for scaler order
				float** S  = new float*[numXBox];
				int**	ns = new int*[numXBox];
				for(int i = 0; i < numXBox; i++){
					S[i]  = new float[numYBox];
					ns[i] = new int[numYBox];
					for(int j = 0; j < numYBox; j++){
						S[i][j]  = 0;
						ns[i][j] = 0;
					}
				}
		
				// Read in X,Y,Z positions for frame
				std::cout << "Reading Frame " << numFrame <<"\t";
				for(int w = 0; w < numWorms; w++)
					for(int p = 0; p < numPerWorm; p++)
						fin >> charTrash >> x[w][p] >> y[w][p] >> floatTrash;
				
				// Dump corner particles at end of file to trash
				for(int t = 0; t < 4; t++)
					fin >> charTrash >> floatTrash >> floatTrash >> floatTrash;	

				// Load S[][] with y/x for all particle pairs
				for(int w = 0; w < numWorms; w++){
					for(int p = 0; p < numPerWorm - 1; p++){
						const float Xcom = (x[w][p] + x[w][p+1])/2.0f;
						const float Ycom = (y[w][p] + y[w][p+1])/2.0f;
						const int	ibox = int(floor(Xcom/boxSize));
						const int	jbox = int(floor(Ycom/boxSize));
						
						double YdivX = (y[w][p] - y[w][p+1])/(x[w][p] - x[w][p+1]); 
						S[ibox][jbox] += YdivX;
						ns[ibox][jbox]++;
					}
				}	
				
				// Get Average Y/X in box
				for(int i = 0; i < numXBox; i++){
					for(int j = 0; j < numYBox; j++){
						if(ns[i][j]!=0) S[i][j] /= float(ns[i][j]);
						else ns[i][j] = -1; // Flag empty boxes							
					}
				}
				
				// Handle Empty boxes
				for(int i = 0; i < numXBox; i++){
					for( int j = 0; j < numYBox; j++){
						if(ns[i][j]==-1){
							int ip1 = (i+1) % numXBox;
							int jp1 = (j+1) % numYBox;
							int im1 = i-1;
							int jm1 = j-1;
							if(i == 0) im1 = numXBox - 1;
							if(j == 0) jm1 = numYBox - 1;
							
							S[i][j] = (S[ip1][j] + S[i][jp1] + S[im1][j] + S[i][jm1])/4;
						}
					}
				}
				
				// Convert average Y/X to angles
				std::cout << "Analysis Completed" << std::endl;
				for(int i = 0; i < numXBox; i++)
					for( int j = 0; j < numYBox; j++)
						S[i][j] = std::atan(S[i][j]);
				
				fout << numBox << std::endl << "comment" << std::endl;
				for(int i = 0; i < numXBox; i++)
					for( int j = 0; j < numYBox; j++)
						fout << "A " << i << " " << j << " " << 0 << " " << std::sin(2*S[i][j])*std::sin(2*S[i][j]) << std::endl;	
				
				// Handle Dynamic Arrays
				for(int dW = 0; dW < numWorms; dW++){
					delete [] x[dW];
					delete [] y[dW];
				}
				
				for(int di = 0; di < numXBox; di++){
					delete [] S[di];
					delete [] ns[di];
				}
					
				delete [] x;
				delete [] y;
				delete [] S;
				delete [] ns;
				numFrame++;
			}
		}
	}
}