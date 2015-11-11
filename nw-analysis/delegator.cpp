#include "delegator.h"

std::string grabsubargs(int& argc, char *argv[], int& i){
	std::string tskarg = "";
	std::string add = "";
	argc = 0;
	while ((add != ") ") || (add != " ) ")){
		tskarg = tskarg + add;
		add = std::string(argv[i++]) + " ";
		argc++;
	}
	return tskarg;
}

// -----------------------------------------------------------------------------------

delegator::delegator(int argc, char *argv[])
{
	//.. process command line arguments
	for (int i = 1; i < argc; ++i){
		std::string arg = argv[i];
		std::string subarg;
		int subargc = 0;

		if (arg == "--follow( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::FOLLOW, subargc, subarg));
		}
		else if (arg == "--gcorr( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::GCORR, subargc, subarg));
		}
		else if (arg == "--gofr( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::GOFR, subargc, subarg));
		}
		else if (arg == "--gtcorr( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::GTCORR, subargc, subarg));
		}
		else if (arg == "--highlight( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::HIGHLIGHT, subargc, subarg));
		}
		else if (arg == "--histo( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::HISTO, subargc, subarg));
		}
		else if (arg == "--persist( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::PERSIST, subargc, subarg));
		}
		else if (arg == "--pofr( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::POFR, subargc, subarg));
		}
		else if (arg == "--pvt( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::PVT, subargc, subarg));
		}
		else if (arg == "--shift2com( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::SHIFT2COM, subargc, subarg));
		}
		else if (arg == "--svden( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::SVDEN, subargc, subarg));
		}
		else if (arg == "--svt( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::SVT, subargc, subarg));
		}
		else if (arg == "--sxyzc( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::SXYZC, subargc, subarg));
		}
		else if (arg == "--varden( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::VARDEN, subargc, subarg));
		}
		else if ((arg == "-h") || (arg == "--help")) 
			this->displayUsage();
		else this->toDo.push_back(new Task(AnalysisTask::UNKNOWN, 0, "unknown"));
	}
}

// -----------------------------------------------------------------------------------

delegator::~delegator()
{
	for (int i = 0; i < toDo.size(); i++)
		toDo.pop_back();
}
