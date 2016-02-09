#pragma once

// ------------------------------------------------
// include all required files storing 
// the calculation functions.
// 2D
#include "follow.h"
#include "gcorr.h"
#include "gofr.h"
#include "gtcorr.h"
#include "highlight.h"
#include "histo.h"
#include "momentum.h"
#include "persistence.h"
#include "pofr.h"
#include "pvt.h"
#include "shift2com.h"
#include "svden.h"
#include "svt.h"
#include "sxyzc.h"
#include "varden.h"

// 3D
#include "thickness.h"
#include "density.h"

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
	( AnalysisTask::GOFR,		&gofr::calculate )
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
	( AnalysisTask::THICKNESS_3D, &thickness_3d::calculate )
	( AnalysisTask::DENSITY_3D, &density_3d::calculate );