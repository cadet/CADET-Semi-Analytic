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

#ifndef VERSIONINFO_HPP_
#define VERSIONINFO_HPP_

namespace casema
{

    //! \brief Returns the version string of the library
    const char* getVersion();

    //! \brief Returns the git commit hash of the source which was used to build the binaries
    const char* getCommitHash();

} // namespace cadet

#endif  // VERSIONINFO_HPP_
