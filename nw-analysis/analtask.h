#pragma once
#include <vector>
#include <string>
#include "boost\thread.hpp"
/*
*	Defines the types of calculations to easy through them
*	around and use.
*/
enum AnalysisTask
{
	FOLLOW,
	GCORR,
	GOFR,
	GTCORR,
	HIGHLIGHT,
	HISTO,
	MOMENTUM,
	PERSIST,
	POFR,
	PVT,
	SHIFT2COM,
	SVDEN,
	SVT,
	SXYZC,
	VARDEN,
	UNKNOWN
};

//.. Groups the thread and type of task together
struct Task
{
	boost::thread thd;
	AnalysisTask tsk;
	std::string arg;
	int argcount;

	Task(AnalysisTask task, int argc = 0, std::string args = ""){
		tsk = task;
		arg = args;
	}
};

//.. define a type for storing tasks
typedef std::vector<Task* > Tasks;