// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_TCLAP_CUSTOM_OUTPUT_HPP_
#define CASEMA_TCLAP_CUSTOM_OUTPUT_HPP_

#include <tclap/StdOutput.h>
#include <string>

#include "VersionInfo.hpp"

namespace TCLAP 
{

	class CustomOutput : public StdOutput
	{
	public:

		CustomOutput(const std::string& progName) : _progName(progName) { }

		virtual void version(CmdLineInterface& c)
		{
			std::cout << "This is " << _progName << " version " << casema::getVersion() << " built from commit " << casema::getCommitHash() << std::endl;
			std::cout << "Report bugs to: cadet@fz-juelich.de" << std::endl;
			std::cout << "CADET homepage: <http://www.cadet-web.de>" << std::endl;
			std::cout << "Fork CADET-semi-analytic on GitHub: <https://github.com/modsim/CADET-semi-analytic>" << std::endl;
		}

	protected:
		std::string _progName;
	};

}

#endif 
