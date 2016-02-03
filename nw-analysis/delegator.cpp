#include "delegator.h"
#include <sstream>
// -----------------------------------------------------------------------------------
//.. include all the analysis functions
#include "analysis_funcs.h"

// -----------------------------------------------------------------------------------
// Returns a string of the char * array starting at index i and updates \p argc, 
// the count of sub strings present until character ')' or string ") " is found.
std::string grabsubargs(int& argc, char *argv[], int& i){
	std::string tskarg = "";
	std::string add = "";
	argc = 0;
	while ((add = std::string(argv[++i])) != "--"){
		tskarg = tskarg + add + " ";
		argc++;
	}
	return tskarg;
}

// -----------------------------------------------------------------------------------
// Returns a vector of all the sub strings in \p str.
std::vector<std::string> parseIntoSubStrings(std::string& str){
	std::vector<std::string> subs;
	std::istringstream iss(str);
	std::string sub;
	while (iss >> sub) subs.push_back(sub);
	return subs;
}

// -----------------------------------------------------------------------------------
// Parses through raw command line arguments and forms container of new tasked to 
// be execuated.  Only access touching parenthesis syntax.
delegator::delegator(int argc, char *argv[]){
	
	//.. print what was recieved from cmdline
	printf("\nCommand line arguments receivced:");
	for (int i = 0; i < argc; i++) printf("\n%i. %s", i, argv[i]);
	
	//.. process command line arguments
	for (int i = 1; i < argc; ++i){
		std::string arg = argv[i];
		std::string subarg;
		int subargc = 0;

		if (arg == "--follow"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::FOLLOW, subargc, subarg));
		}
		else if (arg == "--gcorr"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::GCORR, subargc, subarg));
		}
		else if (arg == "--gofr"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::GOFR, subargc, subarg));
		}
		else if (arg == "--gtcorr"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::GTCORR, subargc, subarg));
		}
		else if (arg == "--highlight"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::HIGHLIGHT, subargc, subarg));
		}
		else if (arg == "--histo"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::HISTO, subargc, subarg));
		}
		else if (arg == "--momentum"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::MOMENTUM, subargc, subarg));
		}
		else if (arg == "--persist"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::PERSIST, subargc, subarg));
		}
		else if (arg == "--pofr"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::POFR, subargc, subarg));
		}
		else if (arg == "--pvt"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::PVT, subargc, subarg));
		}
		else if (arg == "--shift2com"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::SHIFT2COM, subargc, subarg));
		}
		else if (arg == "--svden"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::SVDEN, subargc, subarg));
		}
		else if (arg == "--svt"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::SVT, subargc, subarg));
		}
		else if (arg == "--sxyzc"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::SXYZC, subargc, subarg));
		}
		else if (arg == "--varden"){
			subarg = grabsubargs(subargc, argv, i);
			this->toDo.push_back(new Task(AnalysisTask::VARDEN, subargc, subarg));
		}
		else if ((arg == "-h") || (arg == "--help")) 
			this->displayUsage();
		else this->toDo.push_back(new Task(AnalysisTask::UNKNOWN, 0, "unknown"));
	}
}

// -----------------------------------------------------------------------------------
// Destorys all toDo tasks.
delegator::~delegator(){
	for (int i = 0; i < toDo.size(); i++)
		toDo.pop_back();
}

// -----------------------------------------------------------------------------------
// Assigns all queued tasks to a boost thread, runs them, then waits for them to
// complete.  NOTE:  return values from threads is not supported yet.
int delegator::executeAllTasks(){
	this->assignThreadedFunctions();
	this->joinThreadedFunctions();
	return EXIT_SUCCESS;
}

// -----------------------------------------------------------------------------------
// Traverses through toDo list and initiates the boost::threads in each task with
// the proper calculate function reference and parameters.
void delegator::assignThreadedFunctions(){
	for (auto it : this->toDo){
		if (it->tsk != AnalysisTask::UNKNOWN)
			it->thd = boost::thread(Implement.at(it->tsk), 
							parseIntoSubStrings(it->arg));

		/*switch (it->tsk){
		case AnalysisTask::FOLLOW :
			it->thd = boost::thread(&follow::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << follow::funcName <<"' thread.\n";
			break;
		case AnalysisTask::GCORR :
			it->thd = boost::thread(&gcorr::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << gcorr::funcName << "' thread.\n";
			break;
		case AnalysisTask::GOFR :
			it->thd = boost::thread(&gofr::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << gofr::funcName << "' thread.\n";
			break;
		case AnalysisTask::GTCORR :
			it->thd = boost::thread(&gtcorr::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << gtcorr::funcName << "' thread.\n";
			break;
		case AnalysisTask::HIGHLIGHT :
			it->thd = boost::thread(&highlight::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << highlight::funcName << "' thread.\n";
			break;
		case AnalysisTask::HISTO :
			it->thd = boost::thread(&histo::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << histo::funcName << "' thread.\n";
			break;
		case AnalysisTask::MOMENTUM :
			it->thd = boost::thread(&momentum::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << momentum::funcName << "' thread.\n";
			break;
		case AnalysisTask::PERSIST :
			it->thd = boost::thread(&persistence::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << persistence::funcName << "' thread.\n";
			break;
		case AnalysisTask::POFR :
			it->thd = boost::thread(&pofr::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << pofr::funcName << "' thread.\n";
			break;
		case AnalysisTask::PVT :
			it->thd = boost::thread(&pvt::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << pvt::funcName << "' thread.\n";
			break;
		case AnalysisTask::SHIFT2COM :
			it->thd = boost::thread(&shift2com::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << shift2com::funcName << "' thread.\n";
			break;
		case AnalysisTask::SVDEN :
			it->thd = boost::thread(&svden::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << svden::funcName << "' thread.\n";
			break;
		case AnalysisTask::SVT :
			it->thd = boost::thread(&svt::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << svt::funcName << "' thread.\n";
			break;
		case AnalysisTask::SXYZC :
			it->thd = boost::thread(&sxyzc::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << sxyzc::funcName << "' thread.\n";
			break;
		case AnalysisTask::VARDEN :
			it->thd = boost::thread(&varden::calculate, parseIntoSubStrings(it->arg));
			std::cout << "launching '" << varden::funcName << "' thread.\n";
			break;
		case AnalysisTask::UNKNOWN :
			std::cout << "cannot launch thread for 'unknown' typed task.\n";
			break;
		default :
			std::cout << "default thread launch is not supported.\n";
		}*/
	}
}

// -----------------------------------------------------------------------------------
// Asks all toDo threads to join.
void delegator::joinThreadedFunctions(){
	for (auto it : this->toDo) it->thd.join();
}

// -----------------------------------------------------------------------------------
// Prints to console the usage of this nw-analysis program.
void delegator::displayUsage(){
	std::cout << "Usage:\n\n";
	std::cout << "example:\t>> 'nw-analysis --follow -i in -o out.csv -t 100 -s 1 -e 3 !!'\n\n";
	std::cout << "Commands to the process 'follow' must be followed by '!!' \n";
	std::cout << "e.g. '--process [-option value] ... !!'\n";
}