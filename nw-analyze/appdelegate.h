/*
	This class is used to launch analysis based on user input.

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#pragma once

#ifndef __NW_ANALYSIS_DELEGATE_H__
#define __NW_ANALYSIS_DELEGATE_H__

class nwAnalysisDelegate
{
public:

	//.. only constructor
	nwAnalysisDelegate(int argc, char* argv[]);
	~nwAnalysisDelegate();

	int Run();
private:

	void InterpCLArgs(int argc, char* argv[]);
};

#endif