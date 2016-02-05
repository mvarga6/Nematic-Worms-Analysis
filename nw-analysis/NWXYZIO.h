// Nematic Worms Simulation
// 2.4.16
// Library for versioned i/o to .xyz files
// Mike Varga
#ifndef __NW_UTIL_XYZIO_H__
#define __NW_UTIL_XYZIO_H__
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "NWWormsParameters.h"
#include "NWSimulationParameters.h"
// --------------------------------------
#define isBetween(A, B, C) ( ((A-B) >= -0.0f) && ((A-C) < 0.0f) )
#define VERSION_ERROR_NO -1.0f
// --------------------------------------
// namespace for nematic worms simulation
namespace nw {

	// -----------------------
	// namespace for utilities
	namespace util {

		// ---------------------
		// doing .xyz operations
		namespace xyz {

			// ---------------------------------------------
			// load parameters for version from commnet line
			void readParameterLine(WormsParameters *wparams,
				SimulationParameters *sparams,
				const std::string& comment){

				//---------------------
				//.. get version number
				float version;
				std::string::size_type pos = comment.find_first_of(':');
				if (pos != std::string::npos)
					version = std::strtof(comment.substr(0, pos).c_str(), NULL);
				else version = VERSION_ERROR_NO;

				//---------------------------------
				//.. check for version format error
				if (version == VERSION_ERROR_NO){
					printf("\n*****\nXYZ comment line format error!\n*****");
					return;
				}

				//---------------------------------------
				//.. setup reading stream and tmp storage
				std::stringstream ssvals;
				ssvals << comment.substr(pos, comment.size());
				std::vector<float> fvals;

				// ------------------------------------------------
				//.. decide what to do with version number
				if (isBetween(version, 1.0f, 2.0f)){
					printf("\n%s is undefined for version: %i",
						"nw::util::xyz::readParameterLine",
						(int)version);
				}
				// -------------------------------------------------
				else if (isBetween(version, 2.0f, 3.0f)){

					//.. grab values
					float val;
					while (ssvals >> val) fvals.push_back(val);

					//.. set parameters
					wparams->_NP = (int)fvals.at(0);
					wparams->_Ka = fvals.at(1);
					wparams->_KBT = fvals.at(2);
					wparams->_DRIVE = fvals.at(3);
					wparams->_EPSILON = fvals.at(4);
					wparams->_RCUT = fvals.at(5);
					sparams->_DT = fvals.at(6);
					sparams->_FRAMERATE = (int)fvals.at(7);
				}
				// --------------------------------------------------
				else if (isBetween(version, 3.0f, 4.0f)){

					//.. grab values
					float val;
					while (ssvals >> val) fvals.push_back(val);

					//.. set parameters
					wparams->_NP = (int)fvals.at(0);
					wparams->_Ka = fvals.at(1);
					wparams->_KBT = fvals.at(2);
					wparams->_DRIVE = fvals.at(3);
					wparams->_EPSILON = fvals.at(4);
					wparams->_RCUT = fvals.at(5);
					sparams->_DT = fvals.at(6);
					sparams->_FRAMERATE = (int)fvals.at(7);
				}
				// --------------------------------------------------
				else {
					printf("\n%s is undefined for version: %i",
						"nw::util::xyz::readParameterLine",
						(int)version);
				}

			}

			// -----------------------------------------------
			// construct string of parameters for comment line
			std::string makeParameterLine(WormsParameters *wparams,
				SimulationParameters *sparams,
				const float &version){

				// ----------------------------
				//.. write version number first
				std::stringstream sscom;
				sscom << version << ": ";

				if (isBetween(version, 1.0f, 2.0f)){
					printf("\n%s is undefined for version: %i",
						"nw::util::xyz::makeParameterLine",
						(int)version);
					sscom << "-1.0: " << "\n";
				}
				// -------------------------------------------------
				else if (isBetween(version, 2.0f, 3.0f)){

					//.. write parameters
					sscom << wparams->_NP << ' '
						<< wparams->_Ka << ' '
						<< wparams->_KBT << ' '
						<< wparams->_DRIVE << ' '
						<< wparams->_EPSILON << ' '
						<< wparams->_RCUT << ' '
						<< sparams->_DT << ' '
						<< sparams->_FRAMERATE << "\n";
				}
				// --------------------------------------------------
				else if (isBetween(version, 3.0f, 4.0f)){

					//.. write parameters
					sscom << wparams->_NP << ' '
						<< wparams->_Ka << ' '
						<< wparams->_KBT << ' '
						<< wparams->_DRIVE << ' '
						<< wparams->_EPSILON << ' '
						<< wparams->_RCUT << ' '
						<< sparams->_DT << ' '
						<< sparams->_FRAMERATE << "\n";
				}
				// --------------------------------------------------
				else {
					printf("\n%s is undefined for version: %i",
						"nw::util::xyz::makeParameterLine",
						(int)version);
					sscom << "-1.0: " << "\n";
				}
				return (sscom.str());
			}
		}

	}
}

#endif