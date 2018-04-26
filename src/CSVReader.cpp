// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2018: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "CSVReader.hpp"
#include "StringUtil.hpp"

#include <fstream>

namespace casema
{

	CSVReader::~CSVReader()
	{
		for (std::vector<std::vector<mpfr::mpreal>*>::iterator it = _columnData.begin(); it != _columnData.end(); ++it)
			delete (*it);
	}

	bool CSVReader::contains(const std::string& name) const
	{
		for (std::vector<std::string>::const_iterator it = _columnNames.begin(); it != _columnNames.end(); ++it)
		{
			if (*it == name)
				return true;
		}
		return false;
	}
	
	std::size_t CSVReader::indexOf(const std::string& name) const
	{
		for (std::size_t i = 0; i < _columnNames.size(); ++i)
		{
			if (_columnNames[i] == name)
				return i;
		}
		return -1;
	}

	bool CSVReader::readFile(const std::string& fileName)
	{
		// Open file
		std::ifstream fs(fileName);

		if (!fs.good())
		{
			return false;
		}

		bool header = true;
		std::string line;
		while (std::getline(fs, line))
		{
			util::trim(line);

			// Ignore comments
			if (line[0] == '#')
				continue;

			// Ignore empty lines
			if (line.empty())
				continue;

			// Split line into tokens separated by comma
			const std::vector<std::string> tokens = util::split(line, ",");

			// Read header
			if (header)
			{
				header = false;
				_columnNames = tokens;

				// Allocate vectors for columns
				for (std::size_t i = 0; i < tokens.size(); ++i)
				{
					_columnData.push_back(new std::vector<mpfr::mpreal>());
				}
				continue;
			}

			// Try to convert all items
			for (std::size_t i = 0; i < tokens.size(); ++i)
			{
				const std::size_t prec = mpfr::mpreal::get_default_prec();
				const std::size_t inprec = mpfr::digits2bits(tokens[i].size()+2);

				const char* const str = tokens[i].c_str();
				mpfr::mpreal data(0, std::max(prec, inprec));
				data.setNan();
				char* outPos = 0;
				::mpfr_strtofr(data.mpfr_ptr(), str, &outPos, 10, mpfr::mpreal::get_default_rnd());
				if (outPos == str)
				{
					// Invalid
					_columnData[i]->push_back(mpfr::mpreal().setNan());
				}
				else
				{
					// Valid
					_columnData[i]->push_back(data);
				}
			}
		}

		return true;
	}

}
