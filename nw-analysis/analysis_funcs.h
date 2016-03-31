#pragma once

// ------------------------------------------------
// include all required files storing 
// the calculation functions.
#include "follow.h"			// 2D
#include "gcorr.h"			// 2D
#include "gofr.h"			// 2D & 3D
#include "gtcorr.h"			// 2D
#include "highlight.h"		// 2D
#include "histo.h"			// 2D
#include "momentum.h"		// 2D
#include "persistence.h"	// 2D
#include "pofr.h"			// 2D
#include "pvt.h"			// 2D
#include "shift2com.h"		// 2D
#include "svden.h"			// 2D
#include "svt.h"			// 2D
#include "sxyzc.h"			// 2D & 3D
#include "varden.h"			// 2D
#include "thickness.h"		// 3D
#include "density.h"		// 3D
#include "msd.h"			// 3D
#include "velocity_map.h"	// 2D & quasi-3D

// -------------------------------------------------
// map the AnalysisTasks to functions
// using function pointers.  Const global
// scope maybe can be replaced by an 
// object which links desired functionality
// to implementation.
#include "analtask.h"
#include <map>
#include <vector>
#include <boost/assign.hpp>

using namespace boost::assign;
typedef int(*Calculation)(std::vector<std::string>);

// ------------------------------------------------------------------------
// actual container called "Implement" such that when an analysis function
// is requrested, symantically all thats needed is a loopup in Implement
static const std::map< AnalysisTask, Calculation > Implement = map_list_of
	( AnalysisTask::FOLLOW,		&follow::calculate )
	( AnalysisTask::GCORR,		&gcorr::calculate )
	( AnalysisTask::GOFR,		&gofr::calculate_2d )
	( AnalysisTask::GTCORR,		&gtcorr::calculate )
	( AnalysisTask::HIGHLIGHT,	&highlight::calculate )
	( AnalysisTask::HISTO,		&histo::calculate )
	( AnalysisTask::MOMENTUM,	&momentum::calculate )
	( AnalysisTask::PERSIST,	&persistence::calculate )
	( AnalysisTask::POFR,		&pofr::calculate )
	( AnalysisTask::PVT,		&pvt::calculate )
	( AnalysisTask::SHIFT2COM,	&shift2com::calculate )
	( AnalysisTask::SVDEN,		&svden::calculate )
	( AnalysisTask::SVT,		&svt::calculate )
	( AnalysisTask::SXYZC,		&sxyzc::calculate )
	( AnalysisTask::VARDEN,		&varden::calculate )
	( AnalysisTask::THICKNESS_3D, &thickness::calculate_3d )
	( AnalysisTask::DENSITY_3D, &density::calculate_3d )
	( AnalysisTask::GOFR_3D,	&gofr::calculate_3d )
	( AnalysisTask::SXYZC_3D,	&sxyzc::calculate_3d )
	( AnalysisTask::MSD_3D,		&msd::calculate_3d)
	( AnalysisTask::VGRID,		&velocity::calculate_3d);