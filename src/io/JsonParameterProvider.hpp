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

#include <json.hpp>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <stack>
#include <string>

#include "CompilerSpecific.hpp"
#include "io/ParameterProvider.hpp"

using json = nlohmann::json;

namespace casema
{
namespace io
{

class JsonParameterProvider : public IParameterProvider
{
public:

	JsonParameterProvider(const char* data) : _root(new json(json::parse(data)))
	{
		_opened.push(_root);
	#ifdef CASEMA_DEBUG
		_scopePath = "/";
	#endif
	}

	JsonParameterProvider(const std::string& data) : _root(new json(json::parse(data)))
	{
		_opened.push(_root);
	#ifdef CASEMA_DEBUG
		_scopePath = "/";
	#endif
	}

	JsonParameterProvider(const json& data) : _root(new json(data))
	{
		_opened.push(_root);
	#ifdef CASEMA_DEBUG
		_scopePath = "/";
	#endif
	}

	JsonParameterProvider(const JsonParameterProvider& cpy)
	{
		_root = new json(*cpy._root);
		_opened = cpy._opened;
	#ifdef CASEMA_DEBUG
		_scopePath = cpy._scopePath;
	#endif
	}

	JsonParameterProvider(JsonParameterProvider&& cpy) CASEMA_NOEXCEPT : _root(cpy._root), _opened(std::move(cpy._opened))
	{
		cpy._root = nullptr;
		cpy._opened = std::stack<json*>();
	#ifdef CASEMA_DEBUG
		_scopePath = std::move(cpy._scopePath);
	#endif
	}

	~JsonParameterProvider() CASEMA_NOEXCEPT
	{
		delete _root;
	}

	JsonParameterProvider& operator=(const JsonParameterProvider& cpy)
	{
		delete _root;

		_root = new json(*cpy._root);
		_opened = cpy._opened;

	#ifdef CASEMA_DEBUG
		_scopePath = cpy._scopePath;
	#endif

		return *this;
	}

	JsonParameterProvider& operator=(JsonParameterProvider&& cpy) CASEMA_NOEXCEPT
	{
		delete _root;
		_opened = std::stack<json*>();
		_root = cpy._root;
		cpy._root = nullptr;
		cpy._opened = std::stack<json*>();

		_opened.push(_root);

	#ifdef CASEMA_DEBUG
		_scopePath = std::move(cpy._scopePath);
	#endif

		return *this;
	}

	double getDouble(const std::string& paramName)
	{
		json& p = _opened.top()->at(paramName);
		if (p.is_array() && (p.size() == 1))
			p = p[0];

		return p.get<double>();
	}

	int getInt(const std::string& paramName)
	{
		json& p = _opened.top()->at(paramName);
		if (p.is_array() && (p.size() == 1))
			p = p[0];

		if (p.is_boolean())
			return p.get<bool>();
		else
			return p.get<int>();
	}

	uint64_t getUint64(const std::string& paramName)
	{
		json& p = _opened.top()->at(paramName);
		if (p.is_array() && (p.size() == 1))
			p = p[0];

		return p.get<uint64_t>();
	}

	bool getBool(const std::string& paramName)
	{
		json& p = _opened.top()->at(paramName);
		if (p.is_array() && (p.size() == 1))
			p = p[0];

		if (p.is_number_integer())
			return p.get<int>();
		else
			return p.get<bool>();
	}

	std::string getString(const std::string& paramName)
	{
		json& p = _opened.top()->at(paramName);
		if (p.is_array() && (p.size() == 1))
			p = p[0];

		return p.get<std::string>();
	}

	std::vector<double> getDoubleArray(const std::string& paramName)
	{
		const json& p = _opened.top()->at(paramName);
		if (!p.is_array())
		{
			return std::vector<double>(1, p.get<double>());
		}

		return p.get<std::vector<double>>();
	}

	std::vector<int> getIntArray(const std::string& paramName)
	{
		const json& p = _opened.top()->at(paramName);
		if (p.is_array())
		{
			if (p.size() == 0)
				return std::vector<int>(0);

			if (p[0].is_boolean())
			{
				const std::vector<bool> d = p.template get<std::vector<bool>>();
				std::vector<int> bd(d.size());
				for (std::size_t i = 0; i < d.size(); ++i)
					bd[i] = d[i];

				return bd;
			}

			return p.template get<std::vector<int>>();
		}
		else
		{
			if (p.is_boolean())
			{
				return std::vector<int>(1, p.template get<bool>());
			}

			return std::vector<int>(1, p.template get<int>());
		}
	}

	std::vector<uint64_t> getUint64Array(const std::string& paramName)
	{
		const json& p = _opened.top()->at(paramName);
		if (!p.is_array())
		{
			return std::vector<uint64_t>(1, p.get<uint64_t>());
		}
		return p.get<std::vector<uint64_t>>();
	}

	std::vector<bool> getBoolArray(const std::string& paramName)
	{
		const json& p = _opened.top()->at(paramName);
		if (p.is_array())
		{
			if (p.size() == 0)
				return std::vector<bool>(0);

			if (p[0].is_number_integer())
			{
				const std::vector<int> d = p.template get<std::vector<int>>();
				std::vector<bool> bd(d.size());
				for (std::size_t i = 0; i < d.size(); ++i)
					bd[i] = d[i];

				return bd;
			}

			return p.template get<std::vector<bool>>();
		}
		else
		{
			if (p.is_number_integer())
			{
				return std::vector<bool>(1, p.template get<int>());
			}

			return std::vector<bool>(1, p.template get<bool>());
		}
	}

	std::vector<std::string> getStringArray(const std::string& paramName)
	{
		const json& p = _opened.top()->at(paramName);
		if (!p.is_array())
		{
			return std::vector<std::string>(1, p.get<std::string>());
		}
		return p.get<std::vector<std::string>>();
	}

	bool exists(const std::string& paramName)
	{
		return _opened.top()->find(paramName) != _opened.top()->end();
	}

	bool isArray(const std::string& paramName)
	{
		return _opened.top()->at(paramName).is_array() && (_opened.top()->at(paramName).size() > 1);
	}

	std::size_t numElements(const std::string& paramName)
	{
		return _opened.top()->at(paramName).size();
	}

	void pushScope(const std::string& scope)
	{
		_opened.push(&_opened.top()->at(scope));

	#ifdef CASEMA_DEBUG
		_scopePath += "/" + scope;
	#endif
	}

	void popScope()
	{
		_opened.pop();

	#ifdef CASEMA_DEBUG
		std::size_t lastIdx = std::string::npos;
		if (_scopePath.back() == '/')
			lastIdx = _scopePath.length() - 2;

		const std::size_t idx = _scopePath.find_last_of('/', lastIdx);
		_scopePath.erase(idx);
	#endif
	}

	void addScope(const std::string& scope)
	{
		if (!exists(scope))
		{
			json j;
			j["blubber"] = 0.0;
			j.erase("blubber");
			(*_opened.top())[scope] = j;
		}
	}

	void set(const std::string& paramName, double val)
	{
		(*_opened.top())[paramName] = val;
	}

	void set(const std::string& paramName, int val)
	{
		(*_opened.top())[paramName] = val;
	}

	void set(const std::string& paramName, uint64_t val)
	{
		(*_opened.top())[paramName] = val;
	}

	void set(const std::string& paramName, bool val)
	{
		(*_opened.top())[paramName] = val;
	}

	void set(const std::string& paramName, char const* val)
	{
		(*_opened.top())[paramName] = std::string(val);
	}

	void set(const std::string& paramName, const std::string& val)
	{
		(*_opened.top())[paramName] = val;
	}

	void set(const std::string& paramName, const std::vector<double>& val)
	{
		(*_opened.top())[paramName] = val;
	}

	void set(const std::string& paramName, const std::vector<int>& val)
	{
		(*_opened.top())[paramName] = val;
	}

	void set(const std::string& paramName, const std::vector<uint64_t>& val)
	{
		(*_opened.top())[paramName] = val;
	}

	void set(const std::string& paramName, const std::vector<std::string>& val)
	{
		(*_opened.top())[paramName] = val;
	}

	void remove(const std::string& name)
	{
		(*_opened.top()).erase(name);
	}

	void copy(const std::string& src, const std::string& dest)
	{
		const json j = (*_opened.top())[src];
		(*_opened.top())[dest] = j;
	}

	void toFile(const std::string& fileName) const
	{
	    std::ofstream ofs(fileName, std::ios::out | std::ios::trunc);
	    ofs << _root->dump(4);
	}

	JsonParameterProvider fromFile(const std::string& fileName)
	{
		std::ifstream ifs(fileName);
		json* root = new json();
		ifs >> (*root);

		return JsonParameterProvider(root);
	}

	inline json* data() { return _root; }
	inline json const* data() const { return _root; }

private:
	JsonParameterProvider() : _root(nullptr)
	{
	#ifdef CASEMA_DEBUG
		_scopePath = "/";
	#endif
	}

	JsonParameterProvider(json* data) : _root(data)
	{
		_opened.push(_root);
	#ifdef CASEMA_DEBUG
		_scopePath = "/";
	#endif
	}

	json* _root;
	std::stack<json*> _opened;

#ifdef CASEMA_DEBUG
	std::string _scopePath;
#endif
};

inline std::ostream& operator<<(std::ostream& out, const JsonParameterProvider& jpp)
{
	out << jpp.data()->dump(4);
	return out;
}

} // namespace io
} // namespace casema
