/*
This file defines all the available ways to write a 
comma seperated value (*.csv) to a data file.

M.Varga
Kent State University
Chemical Physics Dept.
2015
*/

#ifndef __CSV_WRITE_HPP__
#define __CSV_WRITE_HPP__

#include <string>
#include <sstream>
#include <vector>
#include <fstream>

using namespace std;

namespace nw
{
	/*
	*	Takes any number of vectors with names and prints them separated by
	*	commas. Prints number of lines equal to max size in list.
	*/
	template <typename dataType>
	void csv_write(string fileNameBase, vector< vector<dataType> >& list, vector<string>& names)
	{
		//.. number of elements in comma seperated lines
		size_t outer_size = list.size();
		if (outer_size == 0) return;

		//.. elements in each column
		vector<size_t> sizes;
		for (int i = 0; i < outer_size; i++)
		{
			sizes.push_back(list[i].size());
		}

		//.. get line max number
		size_t max_size = 0;
		for (int i = 0; i < outer_size; i++)
		{
			if (sizes[i] > max_size) max_size = sizes[i];
		}

		//.. open file
		ofstream fout(fileNameBase + ".csv", ios::out);

		//.. print names of columns to top line
		if (outer_size == names.size())
		{
			fout << "line,";
			for (int i = 0; i < outer_size; i++)
			{
				fout << names[i] << ",";
			}
			fout << endl;
		}

		//.. loop through all lines
		for (int line = 0; line < max_size; line++)
		{
			//.. start each line with line number
			fout << line + 1 << ",";

			//.. print all data for each column
			for (int i = 0; i < outer_size; i++)
			{
				//.. if it exists
				if (line < sizes[i]) fout << list[i][line];
				fout << ",";
			}
			fout << endl;
		}

		//.. close output file
		fout.close();
	}

	/*
	*	Takes any number of vectors without names and prints them separated by
	*	commas. Prints number of lines equal to max size in list.
	*/
	template <typename dataType>
	void csv_write(string fileNameBase, vector< vector<dataType> >& list)
	{
		//.. formulate generic names list
		vector<string> names;
		for (int i = 0; i < list.size(); i++)
		{
			//.. name by number
			ostringstream name;
			name << "data" << i;
			names.push_back(name.str());
		}

		//.. call master function
		csv_write<dataType>(fileNameBase, list, names);
	}
}

#endif
