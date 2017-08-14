// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2017: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "CaSeMaConfig.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <exception>
#include <algorithm>
#include <fstream>

#include "MPReal.hpp"

#include <tclap/CmdLine.h>

#include "TclapCustomOutput.hpp"
#include "VersionInfo.hpp"

#include "CadetEnumeration.hpp"
#include "hdf5/HDF5Reader.hpp"
#include "xml/XMLReader.hpp"
#include "CSVReader.hpp"

struct FileData
{
	std::vector<mpfr::mpreal> time;
	std::vector<mpfr::mpreal> outlet;

	std::size_t nCol;
	std::size_t nPar;

	std::size_t timeOffset;
};


void checkTime(FileData& ref, FileData& cmp)
{
/*
	for (std::size_t i = 0; i < ref.size(); ++i)
	{
		if (ref.time[i] != cmp.time[i+cmp.timeOffset])
	}
*/
}


std::size_t systemSize(const FileData& ds)
{
	return 2 * ds.nPar * ds.nCol + ds.nCol * 2;
}


void compare(std::ostream& fs, std::vector<FileData>& data, std::size_t precision)
{
	// Determine precision
	if (precision <= 1)
	{
		for (std::vector<FileData>::iterator it = data.begin(); it != data.end(); ++it)
		{
			precision = std::max(precision, std::size_t(it->time[0].getPrecision()));
			precision = std::max(precision, std::size_t(it->outlet[0].getPrecision()));
		}
	}

	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	// Set output flags
	std::ios_base::fmtflags curFlags = fs.flags();
	std::streamsize curPrec = fs.precision();

	fs << "syssize,l1,l2,linf,rel1,rel2,relinf\n";
	fs.flags(std::ios::scientific);
	fs.precision(precision);

	// Calculate reference values
	mpfr::mpreal ref1(0);
	mpfr::mpreal ref2(0);
	mpfr::mpreal refinf(0);

	for (std::size_t i = 0; i < data[0].time.size(); ++i)
	{
		ref1 += mpfr::abs(data[0].outlet[i]);
		ref2 += mpfr::sqr(data[0].outlet[i]);
		refinf = mpfr::max(mpfr::abs(data[0].outlet[i]), refinf);
	}
	ref2 = mpfr::sqrt(ref2);

	// Compare dataset times
	for (std::size_t i = 1; i < data.size(); ++i)
	{
		if (data[0].time.size() != data[i].time.size())
		{
			data[i].timeOffset = data[i].time.size() - data[0].time.size();
		}
		checkTime(data[0], data[i]);

		// Calculate norms
		mpfr::mpreal l1(0);
		mpfr::mpreal l2(0);
		mpfr::mpreal linf(0);

		for (std::size_t j = 0; j < data[0].time.size(); ++j)
		{
			mpfr::mpreal diff = mpfr::abs(data[0].outlet[j] - data[i].outlet[j]);
			l1 += diff;
			l2 += mpfr::sqr(diff);
			linf = mpfr::max(diff, linf);
		}
		l2 = mpfr::sqrt(l2);

		fs << systemSize(data[i]) << "," << l1 << "," << l2 << "," << linf << "," << (l1 / ref1) << "," << (l2 / ref2) << "," << (linf / refinf) << "\n";
	}

	// Reset output flags and clear buffers
	fs << std::endl;
	fs.precision(curPrec);
	fs.flags(curFlags);
}


template <class Reader_t>
void readFile(Reader_t& reader, const std::string& fileName, std::vector<FileData>& dest)
{
	reader.openFile(fileName);

	// Check if CADET output or collected data
	if (reader.exists(casema::e2s(casema::GRP_OUT_SOLUTION)))
	{
		// CADET output
		reader.setGroup(casema::e2s(casema::GRP_OUT_SOLUTION));
		const std::vector<double> time = reader.template vector<double>("SOLUTION_COLUMN_OUTLET_COMP_000");
		const std::vector<double> outlet = reader.template vector<double>("SOLUTION_TIMES");

		reader.setGroup(casema::e2s(casema::GRP_IN_DISCRETIZATION));

		FileData data;
		data.nCol = reader.template scalar<int>(casema::e2s(casema::NCOL));
		data.nPar = reader.template scalar<int>(casema::e2s(casema::NPAR));
		data.time.reserve(time.size());
		data.outlet.reserve(outlet.size());
		data.timeOffset = 0;

		for(std::size_t i = 0; i < time.size(); ++i)
		{
			data.time.push_back(mpfr::mpreal(time[i]));
			data.outlet.push_back(mpfr::mpreal(outlet[i]));
		}

		dest.push_back(data);
	}
	else
	{
		// Condensed output
		std::size_t k = 1;

		while (reader.exists("/" + std::to_string(k)))
		{
			reader.setGroup("/" + std::to_string(k));

			const std::vector<double> time = reader.template vector<double>("outlet");
			const std::vector<double> outlet = reader.template vector<double>("time");
			const std::vector<int> sizes = reader.template vector<int>("SYSTEMSIZE");

			FileData data;
			data.time.reserve(time.size());
			data.outlet.reserve(outlet.size());
			data.nCol = sizes[0];
			data.nPar = sizes[1];
			data.timeOffset = 0;

			for(std::size_t i = 0; i < time.size(); ++i)
			{
				data.time.push_back(mpfr::mpreal(time[i]));
				data.outlet.push_back(mpfr::mpreal(outlet[i]));
			}

			dest.push_back(data);

			++k;
		}
	}

	reader.closeFile();
}


int main(int argc, char** argv)
{
	std::size_t precision;
	std::vector<std::string> files;
	std::string outFile;

	try
	{
		TCLAP::CustomOutput customOut("chromError");
		TCLAP::CmdLine cmd("Compares different GRM solutions in arbitrary precision", ' ', casema::getVersion());
		cmd.setOutput(&customOut);

		cmd >> (new TCLAP::ValueArg<std::size_t>("p", "prec", "Precision (default: infer from input, stay below " + std::to_string(MPFR_PREC_MAX) + ")", false, 0, "Int"))->storeIn(&precision);
		cmd >> (new TCLAP::ValueArg<std::string>("o", "out", "Write output to file (default: disabled)", false, std::string(), "File"))->storeIn(&outFile);
		cmd >> (new TCLAP::UnlabeledMultiArg<std::string>("data", "Files", true, "", "File"))->storeIn(&files);

		cmd.parse( argc, argv );
	}
	catch (const TCLAP::ArgException &e)
	{
		std::cerr << "ERROR: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}

	std::vector<FileData> data;
	for (std::vector<std::string>::iterator it = files.begin(); it != files.end(); ++it)
	{
		// Depending on the file type use different reader
		std::size_t fileEnding = it->find_last_of('.') + 1;

		try
		{
			std::string ext = it->substr(fileEnding);
			std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
			if ((ext == "h5") || (ext == "hdf5"))
			{
				casema::HDF5Reader reader;
				readFile(reader, *it, data);
			}
			else if (ext == "xml")
			{
				casema::XMLReader reader;
				readFile(reader, *it, data);
			}
			else if (ext == "csv")
			{
				casema::CSVReader reader;
				if (reader.readFile(*it))
				{
					FileData fd;
					fd.timeOffset = 0;
					fd.nCol = 0;
					fd.nPar = 0;
					fd.time = reader.column(reader.indexOf("time"));
					fd.outlet = reader.column(reader.indexOf("outlet"));
					data.push_back(fd);
				}
				else
					std::cout << "ERROR reading " + (*it) << std::endl;
			}
			else
				std::cout << "ERROR: Wrong input file format \"" + (*it) + "\"! Aborting..." << std::endl;
		}
		catch (const std::exception& e)
		{
			std::cout << "ERROR reading " + (*it) + ": " << e.what() << std::endl;
		}
	}

	std::cout << "Read " << data.size() << " datasets" << std::endl;

	if (data.size() <= 1)
	{
		std::cout << "Please specify more than one dataset" << std::endl;
	}
	else
	{
		if (outFile.empty())
		{
			compare(std::cout, data, precision);
		}
		else
		{
			std::ofstream fs(outFile, std::ofstream::out | std::ofstream::trunc);
			compare(fs, data, precision);
		}
	}

	::mpfr_free_cache();

	return 0;
}

