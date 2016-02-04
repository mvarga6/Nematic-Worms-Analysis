#pragma once
#include <vector>
#include <string>
#include <boost/thread.hpp>

//	Defines the types of calculations to easy through them
//	around and use.
enum AnalysisTask
{
	THICKNESS_3D,
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

// Groups the thread and type of task together
struct Task
{
	//.. thread to execute calculate within
	boost::thread thd;

	//.. type of task requested by user at cmdline
	AnalysisTask tsk;

	//.. arguments for task calculate
	std::string arg;

	//.. count of arguments to task
	int argcount;
	
	//.. init with task type, arg count, and string of args
	Task(AnalysisTask task, int argc = 0, std::string args = ""){
		tsk = task;
		arg = args;
	}
};

// define a type for storing tasks
typedef std::vector<Task* > Tasks;