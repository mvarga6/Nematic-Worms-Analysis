#include "delegator.h"

// COPY COPY COPY

// Create an app delegator with commmand line arguments
// then execute them all on return
int main(int argc, char* argv[]){
	delegator * appDelegate = new delegator(argc, argv);
	return appDelegate->executeAllTasks();
}