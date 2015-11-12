#pragma once

#include <string>
#include <iostream>
#include "analtask.h"
/*
*	Used to interpret cmdline arguments and throw otf Boost threads
*	to run the desired calculations.
*/
class delegator
{
	//.. state variables
	Tasks toDo;

public:
	delegator(int argc, char *argv[]);
	~delegator();

	int executeAllTasks();

private:

	void assignThreadedFunctions();
	void displayUsage();
};

