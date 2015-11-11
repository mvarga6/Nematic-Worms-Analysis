#include "delegator.h"

std::string grabsubargs(char *argv[], int& i){
	std::string tskarg = "";
	std::string add = "";
	while ((add != ") ") || (add != " ) "){
		tskarg = tskarg + add;
		add = std::string(argv[i++]) + " ";
	}
	return tskarg;
}


delegator::delegator(int argc, char *argv[])
{
	//.. process command line arguments
	for (int i = 1; i < argc; ++i){
		std::string arg = argv[i];
		if (arg == "--follow( ") this->toDo.push_back(new Task(AnalysisTask::FOLLOW, grabsubargs(argv, i)));
		else if (arg == "--gcorr( ") this->toDo.push_back(new Task(AnalysisTask::GCORR, grabsubargs(argv, i)));
		else if (arg == "--gofr") this->toDo.push_back(new Task(AnalysisTask::GOFR, grabsubargs(argv, i)));
		else if (arg == "--gtcorr") this->toDo.push_back(new Task(AnalysisTask::GTCORR, grabsubargs(argv, i)));
		else if (arg == "--highlight") this->toDo.push_back(new Task(AnalysisTask::HIGHLIGHT, grabsubargs(argv, i)));
		else if (arg == "--histo") this->toDo.push_back(new Task(AnalysisTask::HISTO, grabsubargs(argv, i)));
		else if (arg == "--persist") this->toDo.push_back(new Task(AnalysisTask::PERSIST, grabsubargs(argv, i)));
		else if (arg == "--pofr") this->toDo.push_back(new Task(AnalysisTask::POFR, grabsubargs(argv, i)));
		else if (arg == "--pvt") this->toDo.push_back(new Task(AnalysisTask::PVT, grabsubargs(argv, i)));
		else if (arg == "--shift2com") this->toDo.push_back(new Task(AnalysisTask::SHIFT2COM, grabsubargs(argv, i)));
		else if (arg == "--svden") this->toDo.push_back(new Task(AnalysisTask::SVDEN, grabsubargs(argv, i)));
		else if (arg == "--svt") this->toDo.push_back(new Task(AnalysisTask::SVT, grabsubargs(argv, i)));
		else if (arg == "--sxyzc") this->toDo.push_back(new Task(AnalysisTask::SXYZC, grabsubargs(argv, i)));
		else if (arg == "--varden") this->toDo.push_back(new Task(AnalysisTask::VARDEN, grabsubargs(argv, i)));
		else if ((arg == "-h") || (arg == "--help")) 
			this->displayUsage();
		else this->toDo.push_back(new Task(AnalysisTask::UNKNOWN, "unknown"));
	}
}


delegator::~delegator()
{
	for (int i = 0; i < toDo.size(); i++)
		toDo.pop_back();
}
