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

#ifndef CASEMA_MODELREADER_HELPER_HPP_
#define CASEMA_MODELREADER_HELPER_HPP_

#include <string>
#include <algorithm>
#include <iostream>

#include "ModelReader.hpp"
#include "hdf5/HDF5Reader.hpp"
#include "xml/XMLReader.hpp"

namespace casema
{
    namespace detail
    {
        template <typename real_t, class reader_t>
        ModelData<real_t> readModel(const std::string& fileName)
        {
            reader_t reader;

            reader.openFile(fileName);
            casema::ModelData<real_t> model = casema::readModel<reader_t, real_t>(reader);
            reader.closeFile();

            return model;
        }
    }

    template <typename real_t>
    bool readModel(const std::string& fileName, ModelData<real_t>& out)
    {
        // Extract suffix
        std::size_t fileEnding = fileName.find_last_of('.') + 1;
        std::string ext = fileName.substr(fileEnding);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        try
        {
            // Select reader
            if (ext == "h5")
                out = detail::readModel<real_t, casema::HDF5Reader>(fileName);
            else if (ext == "xml")
                out = detail::readModel<real_t, casema::XMLReader>(fileName);
            else
                return false;
            
            return true;
        }
        catch (const std::exception& e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
        }

        return false;
    }

}

#endif
