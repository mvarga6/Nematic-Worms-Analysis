
#include "delegator.h"

int main(int argc, char* argv[]){

	delegator * appDelegate = new delegator(argc, argv);
	return appDelegate->executeAllTasks();
}