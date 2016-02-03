#pragma once

#include <boost/thread.hpp>
#include <string>
#include <iostream>
#include "analtask.h"

// Used to interpret cmdline arguments and throw off boost::threads
// to run the desired calculations asyncronously.
class delegator
{
	//.. state variables
	Tasks toDo;

public:
	//.. all that is necessary is to init with cmdline then execute
	delegator(int argc, char *argv[]);
	~delegator();

	//.. execute all tasks queued in constructor
	int executeAllTasks();

private:

	//.. parses through toDo list and throws tasks to threads
	void assignThreadedFunctions();

	//.. asks all tasks to join
	void joinThreadedFunctions();

	//.. print to console this programs usage
	void displayUsage();
};

