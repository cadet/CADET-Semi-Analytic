// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2019: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_CSVREADER_HPP_
#define CASEMA_CSVREADER_HPP_

#include <string>
#include <vector>
#include "MPReal.hpp"

namespace casema
{
	class CSVReader
	{
	public:
		CSVReader() { }
		~CSVReader();

		bool readFile(const std::string& fileName);

		bool contains(const std::string& name) const;
		std::size_t indexOf(const std::string& name) const;

		std::size_t numColumns() const { return _columnData.size(); }

		const std::vector<std::string>& header() const { return _columnNames; }
		const std::string& header(std::size_t idx) const { return _columnNames[idx]; }

		const std::vector<mpfr::mpreal>& column(const std::string& name) const { return *_columnData[indexOf(name)]; }
		const std::vector<mpfr::mpreal>& column(std::size_t idx) const { return *_columnData[idx]; }
		const mpfr::mpreal& last(std::size_t idx) const { return *(_columnData[idx]->end() - 1); }

	protected:
		std::vector<std::string> _columnNames;
		std::vector<std::vector<mpfr::mpreal>*> _columnData;
	};
}

#endif
