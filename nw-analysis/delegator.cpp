#include "delegator.h"
#include <sstream>
#include "follow.h"
#include "gcorr.h"

// -----------------------------------------------------------------------------------
// Returns a string of the char * array starting at index i and updates \p argc, 
// the count of sub strings present until character ')' or string ") " is found.
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
// Returns a vector of all the sub strings in \p str.
std::vector<std::string> parseintosubstrings(std::string& str){
	std::vector<std::string> subs;
	std::istringstream iss(str);
	std::string sub;
	while (iss >> sub) subs.push_back(sub);
	return subs;
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
		else if (arg == "--momentum( "){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::MOMENTUM, subargc, subarg));
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

// -----------------------------------------------------------------------------------

void delegator::assignThreadedFunctions(){
	for (auto it : this->toDo){
		switch (it->tsk){
		case AnalysisTask::FOLLOW :
			it->thd = boost::thread(&follow::calculate, parseintosubstrings(it->arg));
			std::cout << "launching 'follow' thread.\n";
			break;
		case AnalysisTask::GCORR :
			it->thd = boost::thread(&gcorr::calculate, parseintosubstrings(it->arg));
			std::cout << "launching 'gcorr' thread.\n";
		}
	}
}
