#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

//#include "nwAnalysisDelegate.h"

#include "csv_write.hpp"
#include "xyz_read_frame.hpp"
#include "simreplay_types.hpp"
//#include "litterbox.hpp"
#include "simreplay_convert.hpp"
//#include "orientation_hist.hpp"
#include "spatial_correlator.hpp"
//#include "orientation_correlator.hpp"
#include "pairs_data.hpp"
#include "param_set.hpp"
#include "linked_list.hpp"
//#include "xyzv_read_frame.hpp"
//#include "velocity_auto_correlator.hpp"

// DEVELOPMENT MAIN FUNCTION

int main(int argc, char* argv[])
{
	//.. things for startup functions
	// std::srand(time(NULL));
	//  Need to grab first frame before all calculations.
	//	This is setup for things to calculate inter frame
	//	changes, i.e. velocity

	// TEST READING IN AND WRITING 

	std::ifstream fin("SPW1_nwsim_4.xyz");
	nw::parameterSet params("SPW1_cfg.dat");
	std::vector<nw::xyz<float>> xyz;
	//std::vector< std::vector< nw::xyz<float> > > allXYZ;
	nw::xyz_read_frame<float>(fin, xyz);
	//std::vector< std::vector<float>* > allHist;
	//std::vector< std::vector<float>* > allxValues;

	//nw::velocityAutoCorrelator<float> vAutoCorr(params, 20, 0.5f);

	//while (nw::xyz_read_frame<float>(fin, xyz))
	{
		//.. data
		//allXYZ.push_back(xyz);
		//	std::cout << xyz.size() << std::endl;
		//	vAutoCorr.update(allXYZ.back());

		//	vAutoCorr.start();
		//	vAutoCorr.join();

		//.. make histogram
		//std::vector<nw::xyzv<float>> xyzv(nw::xyz_to_xyzv<float>(xyz, params));
		//nw::orientationHistogram<float> hist(xyzv, 45);
		//hist.start();
		//hist.join();

		//.. push onto list of histograms
		//allHist.push_back(new std::vector<float>(hist.results()));
		//allxValues.push_back(new std::vector<float>(hist.xValues()));
	}
	fin.close();

	//for (int i = 0; i < allXYZ.size(); i++)
	//{
	nw::linkedList<float> LL(xyz, params, 10.0f);
	std::vector<nw::xyzv<float>> xyzv(nw::xyz_to_xyzv<float>(xyz, params));
	nw::pairsData<float> pairsInfo(xyzv, 80.0, LL);
	nw::spatialCorrelator<float> spatcorr(xyz, pairsInfo, 2.0f, 60.0f);

	spatcorr.start();
	spatcorr.join();

	//allHist.push_back(new std::vector<float>(spatcorr.results()));
	//allxValues.push_back(new std::vector<float>(spatcorr.xValues()));
	//}


	//.. average
	//std::vector<float> averagedData(nw::averageData<float>(allHist));
	//std::vector<float> averagedXVals(nw::averageData<float>(allxValues));

	//for (int i = 0; i < allHist.size(); i++)
	//{
	//	delete allHist[i];
	//	delete allxValues[i];
	//}

	//nw::linkedList<float> LL(xyz, params, 10.0f);

	//std::vector<nw::xyzv<float>> xyzv(nw::xyz_to_xyzv<float>(xyz, params));

	//nw::orientationHistogram<float> hist(xyzv, 45);
	//hist.start();
	//hist.join();

	std::vector< std::vector<float> > histData = { spatcorr.xValues(), spatcorr.results() };
	//std::vector< std::vector<float> > histData = { averagedXVals, averagedData };
	//std::vector< std::vector<float> > vCorrData = { vAutoCorr.xValues(), vAutoCorr.results() };
	nw::csv_write<float>("spw1_gr", histData);
	//nw::csv_write<float>("AutoCorrelation", vCorrData);

	//nw::pairsData<float> pairsInfo(xyzv, 40.0, LL, 2,false);
	//nw::spatialCorrelator<float> spatcorr(xyz,pairsInfo, 2.0f, 40.0f);
	//nw::orientationCorrelator<float> oriecorr(xyzv, pairsInfo, 2.0f, 40.0f);

	//spatcorr.start();
	//oriecorr.start();
	//spatcorr.join();
	//oriecorr.join();

	//std::vector < std::vector<float> > corrData1 = { spatcorr.xValues(), spatcorr.results() };
	//std::vector < std::vector<float> > corrData2 = { oriecorr.xValues(), oriecorr.results() };

	//nw::csv_write<float>("gr", corrData1);
	//nw::csv_write<float>("nocf", corrData2);

	// BOX TESTING
	//nw::litterBox * box = new nw::periodicBox(10, 12, 14);

	//double x = 5.7, y = 13.2, z = -2.0;
	//double dx = 2.0, dy = 10.0, dz = -8.0;
	//nw::xyz<double> pos = { x, y, z };

	//std::cout << "Before Boundary Conditions" << std::endl;
	//std::cout << pos.x << " " << pos.y << " " << pos.z << std::endl;
	//std::cout << dx << " " << dy << " " << dz << std::endl;

	//box->executeBC(pos);
	//box->recalcDistance(dx, dy, dz);

	//std::cout << "After Boundary Conditions" << std::endl;
	//std::cout << pos.x << " " << pos.y << " " << pos.z << std::endl;
	//std::cout << dx << " " << dy << " " << dz << std::endl;

	// TESTING Convertion xyz -> xyzv;
	//std::ifstream fin("42long_earlier100.xyz");
	//nw::litterBox<float> * box = new nw::periodicBox<float>(10.0, 12.0, 14.0);
	//std::vector<nw::xyz<float>> stuff;
	//nw::xyz_read_frame(fin, stuff);
	//std::vector<nw::xyzv<float>> translated(nw::xyz_to_xyzv(stuff, 42, box));

	//for (auto it : translated)
	//{
	//std::cout << it.t << "\t" << it.x << " " << it.y << " " << it.z
	//<< std::endl << "\t" << it.vx << " " << it.vy << " " << it.vz << std::endl;
	//}


	return 0;
}

// RELEASE ENTRY POINT
/*
int main(int argc, char *argv[])
{
nwAnalysisDelegate * runApp = new nwAnalysisDelegate(argc, argv);
return runApp->Run;
}*/