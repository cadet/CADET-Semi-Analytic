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

#ifndef CASEMA_STRUTIL_HPP_
#define CASEMA_STRUTIL_HPP_

#include <string>
#include <sstream>
#include <vector>
#include <functional> 
#include <cctype>
#include <locale>
#include <algorithm>

namespace casema
{
	namespace util
	{

		inline bool isOneOf(char c, const std::string& delims)
		{
			for (std::size_t i = 0; i < delims.size(); ++i)
			{
				if (c == delims[i])
					return true;
			}
			return false;
		}

		inline std::vector<std::string> split(const std::string& s, const std::string& delims)
		{
			std::vector<std::string> tokens;
			int idxStart = 0;
			bool inToken = false;
			
			for (std::size_t i = 0; i < s.length(); ++i)
			{
				if (isOneOf(s[i], delims))
				{
					if (!inToken)
					{
						idxStart = i+1;
					}
					else
					{
						tokens.push_back(s.substr(idxStart, i - idxStart));
						inToken = false;
						idxStart = i+1;
					}
				}
				else
				{
					inToken = true;
				}
			}

			if (inToken)
			{
				tokens.push_back(s.substr(idxStart));
			}
			return tokens;
		}

		inline std::vector<std::string>& split(const std::string& s, char delim, std::vector<std::string>& elems)
		{
			std::stringstream ss(s);
			std::string item;
			while (std::getline(ss, item, delim))
			{
				elems.push_back(item);
			}
			return elems;
		}

		inline std::vector<std::string> split(const std::string& s, char delim)
		{
			std::vector<std::string> elems;
			split(s, delim, elems);
			return elems;
		}

		inline std::string& ltrim(std::string& s) 
		{
			s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			return s;
		}

		inline std::string& rtrim(std::string& s) 
		{
			s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
			return s;
		}

		inline std::string& trim(std::string& s) 
		{
			return ltrim(rtrim(s));
		}

		inline std::string& toUpper(std::string& s)
		{
			std::transform(s.begin(), s.end(), s.begin(), ::toupper);
			return s;
		}

		inline std::string& toLower(std::string& s)
		{
			std::transform(s.begin(), s.end(), s.begin(), ::tolower);
			return s;
		}

		/**
		 * @brief Checks whether two strings are equal regardless of upper and lower case
		 * @param [in] a First string
		 * @param [in] b Second string
		 * @return @c true if the strings are equal, otherwise @c false
		 */
		inline bool caseInsensitiveEquals(const std::string& a, const std::string& b)
		{
			const unsigned int sz = a.size();
			if (b.size() != sz)
				return false;
			for (unsigned int i = 0; i < sz; ++i)
			{
				if (tolower(a[i]) != tolower(b[i]))
					return false;
			}
			return true;
		}
	}
}

#endif
