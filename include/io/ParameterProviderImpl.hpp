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

/**
 * @file 
 * Provides an implementation of the casema::IParameterProvider interface
 */

#ifndef CASEMA_PARAMPROVIDERIMPL_HPP_
#define CASEMA_PARAMPROVIDERIMPL_HPP_

#include <string>

#include "io/ParameterProvider.hpp"
#include "io/FileIO.hpp"

namespace casema
{
namespace io
{

class FileParameterProvider : public IParameterProvider
{
public:

	FileParameterProvider(IFileReader& reader) : FileParameterProvider(reader, true) { }

	FileParameterProvider(IFileReader& reader, bool inputPrefix) : _reader(reader)
	{
		if (inputPrefix)
			_reader.setGroup("input");
	}

	virtual ~FileParameterProvider() CASEMA_NOEXCEPT { }

	virtual double getDouble(const std::string& paramName)
	{
		return _reader.getDouble(paramName);
	}

	virtual int getInt(const std::string& paramName)
	{
		return _reader.getInt(paramName);
	}

	virtual uint64_t getUint64(const std::string& paramName)
	{
		return _reader.getUint64(paramName);
	}

	virtual bool getBool(const std::string& paramName)
	{
		return _reader.getBool(paramName);
	}

	virtual std::string getString(const std::string& paramName)
	{
		return _reader.getString(paramName);
	}

	virtual std::vector<double> getDoubleArray(const std::string& paramName)
	{
		return _reader.getDoubleArray(paramName);
	}

	virtual std::vector<int> getIntArray(const std::string& paramName)
	{
		return _reader.getIntArray(paramName);
	}

	virtual std::vector<uint64_t> getUint64Array(const std::string& paramName)
	{
		return _reader.getUint64Array(paramName);
	}

	virtual std::vector<bool> getBoolArray(const std::string& paramName)
	{
		return _reader.getBoolArray(paramName);
	}

	virtual std::vector<std::string> getStringArray(const std::string& paramName)
	{
		return _reader.getStringArray(paramName);
	}

	virtual bool exists(const std::string& paramName)
	{
		return _reader.exists(paramName);
	}

	virtual bool isArray(const std::string& paramName)
	{
		return _reader.isArray(paramName);
	}

	virtual std::size_t numElements(const std::string& paramName)
	{
		return _reader.numElements(paramName);
	}

	virtual void pushScope(const std::string& scope)
	{
		_reader.pushGroup(scope);
	}

	virtual void popScope()
	{
		_reader.popGroup();
	}

private:
	IFileReader& _reader;
};

} // namespace io
} // namespace casema

#endif  // CASEMA_PARAMPROVIDERIMPL_HPP_
