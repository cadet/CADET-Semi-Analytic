// =============================================================================
//  CADET-semi-analytic - The semi-analytic extension of CADET
//  
//  Copyright © 2015-2020: Samuel Leweke¹²
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//    ² University of Cologne, Cologne, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef LIBCASEMA_TCLAPUTILS_HPP_
#define LIBCASEMA_TCLAPUTILS_HPP_

#include <tclap/StdOutput.h>
#include <string>
#include <mpreal.h>

#include "VersionInfo.hpp"

namespace TCLAP 
{

	/**
	 * @brief Modifies the standard behavior of TCLAP to output a better version notice
	 * @details The version notice includes the version of CASEMA (tag, branch, and commit hash)
	 *          as well hints to CASEMA-Web and the GitHub project.
	 */
	class CustomOutput : public StdOutput
	{
	public:

		CustomOutput(const std::string& progName) : _progName(progName) { }

		virtual void version(CmdLineInterface& c)
		{
			std::cout << "This is " << _progName << " version " << casema::getVersion() << " (" << casema::getBranchRefspec() << " branch)\n";
			std::cout << "Built from commit " << casema::getCommitHash() << "\n";
		    std::cout << "CASEMA homepage: <http://www.cadet-web.de>\n";
		    std::cout << "Fork CASEMA on GitHub: <https://github.com/modsim/CASEMA>\n";
		    std::cout << "Report bugs to the issue tracker on GitHub or <cadet@fz-juelich.de>" << std::endl;
		}

	protected:
		std::string _progName;
	};


	/**
	 * @brief Modifies the standard behavior of TCLAP to output a better version notice (without library version)
	 * @details The version notice includes hints to CASEMA-Web and the GitHub project.
	 */
	class CustomOutputWithoutVersion : public StdOutput
	{
	public:

		CustomOutputWithoutVersion(const std::string& progName) : _progName(progName) { }

		virtual void version(CmdLineInterface& c)
		{
			std::cout << "This is " << _progName << "\n";
		    std::cout << "CASEMA homepage: <http://www.cadet-web.de>\n";
		    std::cout << "Fork CASEMA on GitHub: <https://github.com/modsim/CASEMA>\n";
		    std::cout << "Report bugs to the issue tracker on GitHub or cadet@fz-juelich.de" << std::endl;
		}

	protected:
		std::string _progName;
	};

	// Enable TCLAP support of mpreal datatypes
	template<>
	struct ArgTraits<mpfr::mpreal>
	{
		typedef StringLike ValueCategory;
	};
	
	template<>
	void SetString(mpfr::mpreal &dst, const std::string &src)
	{
		const std::size_t prec = mpfr::mpreal::get_default_prec();
		const std::size_t inprec = mpfr::digits2bits(src.size()+2);
		
		mpfr::mpreal x(0, std::max(prec, inprec));
		mpfr_strtofr (x.mpfr_ptr(), src.c_str(), 0, 10, mpfr::mpreal::get_default_rnd());

		dst = x;
	}

}

#endif  // LIBCASEMA_TCLAPUTILS_HPP_
