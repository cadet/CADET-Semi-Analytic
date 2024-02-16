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

#ifndef VERSIONINFO_HPP_
#define VERSIONINFO_HPP_

#include "casemaCompilerInfo.hpp"

namespace casema
{

	const char* getVersion() CASEMA_NOEXCEPT;
	const char* getCommitHash() CASEMA_NOEXCEPT;
	const char* getBranchRefspec() CASEMA_NOEXCEPT;
	const char* getDependencyVersions() CASEMA_NOEXCEPT;
	const char* getBuildType() CASEMA_NOEXCEPT;
	const char* getCompiler() CASEMA_NOEXCEPT;
	const char* getCompilerFlags() CASEMA_NOEXCEPT;
	const char* getBuildHost() CASEMA_NOEXCEPT;

} // namespace casema

#endif  // VERSIONINFO_HPP_
