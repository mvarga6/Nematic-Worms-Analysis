#include "param_set.hpp"
#include <sstream>
#include <iostream>

nw::parameterSet::parameterSet()
{
}

nw::parameterSet::parameterSet(std::string fileName)
{
	this->LoadWithName(fileName);
}

nw::parameterSet::~parameterSet()
{
}

bool nw::parameterSet::LoadWithStream(std::ifstream* fin)
{
	//.. stop if invalid stream
	if (!fin->is_open()) return false;

	std::string line;

	//.. get parameters
	while (std::getline(*fin, line))
	{
		if (line == "...DEBUGING...") break;
		std::istringstream linestream(line);
		std::string key;
		if (std::getline(linestream, key, '='))
		{
			std::string value;
			if (std::getline(linestream, value))
			{

				//std::cout << key << "=" << value << std::endl;

				if (key == "version")				set[key].f = std::stof(value);
				else if (key == "key")				simkey = value;
				else if (key == "np")				set[key].i = std::stoi(value);
				else if (key == "nworms")			set[key].i = std::stoi(value);
				else if (key == "nsteps")			set[key].i = std::stoi(value);
				else if (key == "nparticles")		set[key].i = std::stoi(value);
				else if (key == "dt")				set[key].f = std::stof(value);
				else if (key == "mu")				set[key].f = std::stof(value);
				else if (key == "xbox")				set[key].f = std::stof(value);
				else if (key == "ybox")				set[key].f = std::stof(value);
				else if (key == "blocksperkernel")	set[key].i = std::stoi(value);
				else if (key == "threadsperblock")	set[key].i = std::stoi(value);
				else if (key == "epsilon")			set[key].f = std::stof(value);
				else if (key == "sigma")			set[key].f = std::stof(value);
				else if (key == "drive")			set[key].f = std::stof(value);
				else if (key == "criticaldot")		set[key].f = std::stof(value);
				else if (key == "k1")				set[key].f = std::stof(value);
				else if (key == "k2")				set[key].f = std::stof(value);
				else if (key == "k3")				set[key].f = std::stof(value);
				else if (key == "l1")				set[key].f = std::stof(value);
				else if (key == "l2")				set[key].f = std::stof(value);
				else if (key == "l3")				set[key].f = std::stof(value);
				else if (key == "rcut")				set[key].f = std::stof(value);
				else if (key == "framerate")		set[key].i = std::stoi(value);
				else if (key == "nlistmax")			set[key].i = std::stoi(value);
				else if (key == "nlistsetgap")		set[key].i = std::stoi(value);
				else if (key == "extensile")		set[key].i = std::stoi(value);
				else if (key == "noise")			set[key].f = std::stof(value);
				else if (key == "drag")				set[key].f = std::stof(value);
				else if (key == "damp")				set[key].f = std::stof(value);
				else if (key == "sticky")			set[key].i = std::stoi(value);
				else if (key == "fliprate")			set[key].f = std::stof(value);
				else if (key == "flippop")			set[key].f = std::stof(value);
			}
			else
			{
				std::cerr << "No value for '" << key <<"', assuming float = 0.0f ." << std::endl;
				set[key].f = 0.0f;
			}
		}
	}

	return true;
}

bool nw::parameterSet::LoadWithName(std::string fileName)
{
	//.. create stream
	std::ifstream * instream = new std::ifstream(fileName, std::ios::in);

	//.. stop if no file
	if (!instream->is_open()) return false;

	//.. attempt file reading
	bool loadState = LoadWithStream(instream);
	instream->close();

	return loadState;
}

void nw::parameterSet::Show(void)
{
	std::map<std::string, nw::parameter>::iterator it = this->set.begin();
	while (it != this->set.end())
	{
		std::cout << it->first << " = " << it->second.f << std::endl;
		it++;
	}
}

nw::parameter nw::parameterSet::operator[](std::string key) const
{
	return this->set.at(key);
}