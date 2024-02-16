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

#include <catch.hpp>

#include "io/JsonParameterProvider.hpp"
#include "ModelBuilder.hpp"
#include "model/ModelSystem.hpp"
#include "model/UnitOperation.hpp"
#include "ExpInt.hpp"

using json = nlohmann::json;

namespace
{
	json createGeneralRateModel(double colDispersion = 5.75e-8, double colLength = 0.014, double kA = 35.5)
	{
		json j = json::parse(R"json({
			"UNIT_TYPE": "GENERAL_RATE_MODEL",
			"NCOMP": 1,
			"VELOCITY": 5.75e-4,
			"COL_DISPERSION": 5.75e-8,
			"FILM_DIFFUSION": [6.9e-6],
			"PAR_DIFFUSION": [6.07e-11],
			"PAR_SURFDIFFUSION": [0.0],
			"COL_LENGTH": 0.014,
			"PAR_RADIUS": 4.5e-5,
			"COL_POROSITY": 0.37,
			"PAR_POROSITY": 0.75,
			"INIT_C": [0.0],
			"INIT_Q": [0.0],
			"ADSORPTION_MODEL": "LINEAR",
			"adsorption":
			{
				"IS_KINETIC": 0,
				"LIN_KA": [35.5],
				"LIN_KD": [1000.0]
			},
			"discretization":
			{
				"NCOL": 16,
				"NPAR": 4,
				"NBOUND": [1],
				"PAR_DISC_TYPE": "EQUIDISTANT_PAR",
				"USE_ANALYTIC_JACOBIAN": true,
				"MAX_KRYLOV": 0,
				"GS_TYPE": 1,
				"MAX_RESTARTS": 10,
				"SCHUR_SAFETY": 1e-8,
				"weno":
				{
					"WENO_ORDER": 3,
					"BOUNDARY_MODEL": 0,
					"WENO_EPS": 1e-10
				}
			}
			})json");
		j["COL_DISPERSION"] = colDispersion;
		j["COL_LENGTH"] = colLength;
		j["adsorption"]["LIN_KA"] = std::vector<double>{kA};
		return j;
	}

	json createCSTR(double V0)
	{
		json j = json::parse(R"json({
			"UNIT_TYPE": "CSTR",
			"NCOMP": 1,
			"INIT_C": 0.0,
			"INIT_VOLUME": 0.0
			})json");
		j["INIT_VOLUME"] = V0;
		return j;
	}

	json createInlet(double cin)
	{
		json j = json::parse(R"json({
			"UNIT_TYPE": "INLET",
			"NCOMP": 1,
			"sec_000":
			{
				"CONST_COEFF": [0],
				"LIN_COEFF": [0],
				"QUAD_COEFF": [0],
				"CUBE_COEFF": [0]
			},
			"sec_001":
			{
				"CONST_COEFF": [0],
				"LIN_COEFF": [0],
				"QUAD_COEFF": [0],
				"CUBE_COEFF": [0]
			}
			})json");
		j["sec_000"]["CONST_COEFF"] = std::vector<double>{cin};
		return j;
	}

	json createOutlet()
	{
		return json::parse(R"json({
			"UNIT_TYPE": "OUTLET",
			"NCOMP": 1
			})json");
	}

	mpfr::mpreal constantFromUnit(const json& jpp, const mpfr::mpreal& T, const mpfr::mpreal& qIn, const mpfr::mpreal& abscissa, const mpfr::mpreal& injTime)
	{
		const std::string type = jpp["UNIT_TYPE"].get<std::string>();
		if (type == "CSTR")
		{
			return qIn / jpp["INIT_VOLUME"].get<double>() * T / mpfr::const_pi();
		}
		else if (type == "INLET")
		{
			const json& p = jpp["sec_000"];
			const std::vector<double> v = p["CONST_COEFF"].get<std::vector<double>>();
			const mpfr::mpreal cin = v[0];

			return cin * (1 + exp(-abscissa * injTime)) * T / mpfr::const_pi();
		}
		else
		{
			const mpfr::mpreal Dax = jpp["COL_DISPERSION"].get<double>();
			const mpfr::mpreal L = jpp["COL_LENGTH"].get<double>();

			mpfr::mpreal A(-1);
			if (jpp.find("CROSS_SECTION_AREA") != jpp.end())
				A = jpp["CROSS_SECTION_AREA"].get<double>();
			else if (jpp.find("COL_RADIUS") != jpp.end())
				A = mpfr::const_pi() * sqr(mpfr::mpreal(jpp["COL_RADIUS"].get<double>()));

			mpfr::mpreal halfPe;
			if (A > 0.0)
			{
				mpfr::mpreal eps(1);

				if (jpp.find("COL_POROSITY") != jpp.end())
					eps = jpp["COL_POROSITY"].get<double>();
				else if (jpp.find("TOTAL_POROSITY") != jpp.end())
					eps = jpp["TOTAL_POROSITY"].get<double>();

				halfPe = L * qIn / (A * eps) / (2 * Dax);
			}
			else
				halfPe = L * mpfr::mpreal(jpp["VELOCITY"].get<double>()) / (2 * Dax);

			return 8 * sqrt(mpfr::mpreal(2)) * exp(halfPe) * Dax / sqr(exp(mpfr::mpreal(1)) * L) * T / mpfr::const_pi();
		}
	}

	mpfr::mpreal truncationErrorFromUnit(const json& jpp, const mpfr::mpreal& T, const mpfr::mpreal& qIn, const mpfr::mpreal& abscissa, const mpfr::mpreal& M, const mpfr::mpreal& nTerms)
	{
		const std::string type = jpp["UNIT_TYPE"].get<std::string>();
		if (type == "CSTR")
		{
			const mpfr::mpreal V0 = jpp["INIT_VOLUME"].get<double>();
			const mpfr::mpreal qOutAV0 = qIn + V0 * abscissa;
			return 2 * M * T * qIn / (mpfr::const_pi() * qOutAV0) * asinh(qOutAV0 * T / (mpfr::const_pi() * V0 * nTerms));
		}
		else
		{
			const mpfr::mpreal Dax = jpp["COL_DISPERSION"].get<double>();
			const mpfr::mpreal L = jpp["COL_LENGTH"].get<double>();

			mpfr::mpreal A(-1);
			if (jpp.find("CROSS_SECTION_AREA") != jpp.end())
				A = jpp["CROSS_SECTION_AREA"].get<double>();
			else if (jpp.find("COL_RADIUS") != jpp.end())
				A = mpfr::const_pi() * sqr(mpfr::mpreal(jpp["COL_RADIUS"].get<double>()));

			mpfr::mpreal halfPe;
			if (A > 0.0)
			{
				mpfr::mpreal eps(1);

				if (jpp.find("COL_POROSITY") != jpp.end())
					eps = jpp["COL_POROSITY"].get<double>();
				else if (jpp.find("TOTAL_POROSITY") != jpp.end())
					eps = jpp["TOTAL_POROSITY"].get<double>();

				halfPe = L * qIn / (A * eps) / (2 * Dax);
			}
			else
				halfPe = L * mpfr::mpreal(jpp["VELOCITY"].get<double>()) / (2 * Dax);

			const mpfr::mpreal gamma = L * sqrt(mpfr::const_pi() / (2 * T * Dax));
			return sqrt(mpfr::mpreal("32")) * M * T / mpfr::const_pi() * exp(halfPe) * casema::expInt(gamma * sqrt(nTerms));
		}
	}

	json makeConnections(const std::vector<double>& c)
	{
		json conSwitch0;
		conSwitch0["CONNECTIONS"] = c;
		conSwitch0["SECTION"] = 0;

		json con;
		con["NSWITCHES"] = 1;
		con["switch_000"] = conSwitch0;

		return con;
	}
}

TEST_CASE("Time domain upper bound linear chain GRM", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(20));

	const double cIn = 2.0;

	json jsonModel;
	jsonModel["NUNITS"] = 2;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createGeneralRateModel();
	jsonModel["connections"] = makeConnections({ 0.0, 1.0, -1.0, -1.0, 1.0 });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());
	CHECK(sys->timeDomainUpperBound() == cIn);
}

TEST_CASE("Time domain upper bound linear chain 2xGRM", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(20));

	const double cIn = 2.0;

	json jsonModel;
	jsonModel["NUNITS"] = 4;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createGeneralRateModel();
	jsonModel["unit_002"] = createGeneralRateModel();
	jsonModel["unit_003"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 1.0, -1.0, -1.0, 1.0,
	                                             1.0, 2.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, -1.0, -1.0, 1.0});
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());
	CHECK(sys->timeDomainUpperBound() == cIn);
}

TEST_CASE("Time domain upper bound DAG GRM", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(20));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 7;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);
	jsonModel["unit_002"] = createGeneralRateModel();
	jsonModel["unit_003"] = createGeneralRateModel();
	jsonModel["unit_004"] = createGeneralRateModel();
	jsonModel["unit_005"] = createOutlet();
	jsonModel["unit_006"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, -1.0, -1.0, 1.0,
	                                             1.0, 3.0, -1.0, -1.0, 3.0,
	                                             3.0, 4.0, -1.0, -1.0, 4.0,
	                                             4.0, 5.0, -1.0, -1.0, 2.0,
	                                             4.0, 6.0, -1.0, -1.0, 2.0
	                                             });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());
	CHECK(abs(sys->timeDomainUpperBound() - std::max(cIn, cIn2)) <= mpfr::mpreal("1e-10"));
}

TEST_CASE("Time domain upper bound DAG CSTR", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(20));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 7;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);
	jsonModel["unit_002"] = createCSTR(0.1);
	jsonModel["unit_003"] = createCSTR(1.0);
	jsonModel["unit_004"] = createGeneralRateModel();
	jsonModel["unit_005"] = createOutlet();
	jsonModel["unit_006"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, -1.0, -1.0, 1.0,
	                                             1.0, 3.0, -1.0, -1.0, 3.0,
	                                             3.0, 4.0, -1.0, -1.0, 4.0,
	                                             4.0, 5.0, -1.0, -1.0, 2.0,
	                                             4.0, 6.0, -1.0, -1.0, 2.0
	                                           });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());
	CHECK(abs(sys->timeDomainUpperBound() - std::max(cIn, cIn2)) <= mpfr::mpreal("1e-10"));
}

TEST_CASE("Estimate vs inverse GRM", "[Estimate][Unit]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(40));

	const json jsonUnit = createGeneralRateModel();
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonUnit);
	casema::ModelBuilder mb;
	casema::model::UnitOperation* const m = mb.createUnitOperation(jpp, 0);
	m->setSectionTimes(secTimes.data(), secTimes.size());

	std::vector<mpfr::mpreal> flowRates{1};
	m->setFlowRates(flowRates.data(), flowRates.data());

	REQUIRE(m->hasValidEstimate());

	const mpfr::mpreal abscissa(2);
	const mpfr::mpreal T(3);
	const mpfr::mpreal M(7);

	{
		const mpfr::mpreal err("1e-20");
		const mpfr::mpreal nTerms = m->inverseTruncationError(M, abscissa, T, err);
		const mpfr::mpreal truncErr = m->truncationError(M, abscissa, T, nTerms);
		CAPTURE(nTerms);
		CAPTURE(truncErr);
		CHECK(truncErr <= err);
	}

	{
		const mpfr::mpreal nTerms(100);
		const mpfr::mpreal truncErr = m->truncationError(M, abscissa, T, nTerms);
		const mpfr::mpreal nTerms2 = m->inverseTruncationError(M, abscissa, T, truncErr);
		CAPTURE(nTerms2);
		CAPTURE(truncErr);
		CHECK(nTerms2 >= nTerms);
	}
}

TEST_CASE("Estimate vs inverse CSTR", "[Estimate][Unit]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(40));

	const json jsonUnit = createCSTR(1.0);
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonUnit);
	casema::ModelBuilder mb;
	casema::model::UnitOperation* const m = mb.createUnitOperation(jpp, 0);
	m->setSectionTimes(secTimes.data(), secTimes.size());

	std::vector<mpfr::mpreal> flowRates{1};
	m->setFlowRates(flowRates.data(), flowRates.data());

	REQUIRE(m->hasValidEstimate());

	const mpfr::mpreal abscissa(2);
	const mpfr::mpreal T(3);
	const mpfr::mpreal M(7);

	{
		const mpfr::mpreal err("1e-20");
		const mpfr::mpreal nTerms = m->inverseTruncationError(M, abscissa, T, err);
		const mpfr::mpreal truncErr = m->truncationError(M, abscissa, T, nTerms);
		CAPTURE(nTerms);
		CAPTURE(truncErr);
		CHECK(truncErr <= err);
	}

	{
		const mpfr::mpreal nTerms(100);
		const mpfr::mpreal truncErr = m->truncationError(M, abscissa, T, nTerms);
		const mpfr::mpreal nTerms2 = m->inverseTruncationError(M, abscissa, T, truncErr);
		CAPTURE(nTerms2);
		CAPTURE(truncErr);
		CHECK(nTerms2 >= nTerms);
	}
}

TEST_CASE("Truncation error DAG GRM", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(20));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 7;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);
	jsonModel["unit_002"] = createGeneralRateModel();
	jsonModel["unit_003"] = createGeneralRateModel();
	jsonModel["unit_004"] = createGeneralRateModel();
	jsonModel["unit_005"] = createOutlet();
	jsonModel["unit_006"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, -1.0, -1.0, 1.0,
	                                             1.0, 3.0, -1.0, -1.0, 3.0,
	                                             3.0, 4.0, -1.0, -1.0, 4.0,
	                                             4.0, 5.0, -1.0, -1.0, 2.0,
	                                             4.0, 6.0, -1.0, -1.0, 2.0
	                                           });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());

	const mpfr::mpreal abscissa(2);
	const mpfr::mpreal T(3);
	const mpfr::mpreal nTerms(97);

	mpfr::mpreal ref(0);
	{
		const mpfr::mpreal M = (constantFromUnit(jsonModel["unit_000"], T, 0, abscissa, 1) * constantFromUnit(jsonModel["unit_002"], T, 1, abscissa, 1) * 0.25 + constantFromUnit(jsonModel["unit_001"], T, 0, abscissa, 1) * 0.75) * constantFromUnit(jsonModel["unit_003"], T, 4, abscissa, 1);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_004"], T, 4, abscissa, M, nTerms));
	}
	{
		const mpfr::mpreal M = (constantFromUnit(jsonModel["unit_000"], T, 0, abscissa, 1) * constantFromUnit(jsonModel["unit_002"], T, 1, abscissa, 1) * 0.25 + constantFromUnit(jsonModel["unit_001"], T, 0, abscissa, 1) * 0.75);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_003"], T, 4, abscissa, M, nTerms));
	}
	{
		const mpfr::mpreal M = constantFromUnit(jsonModel["unit_000"], T, 0, abscissa, 1);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_002"], T, 1, abscissa, M, nTerms));
	}

	const mpfr::mpreal truncErr = sys->truncationError(abscissa, T, nTerms);
	CHECK(abs((truncErr - ref) / ref) <= mpfr::mpreal("1e-16"));
}

TEST_CASE("Truncation error DAG CSTR", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(20));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 7;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);
	jsonModel["unit_002"] = createCSTR(0.1);
	jsonModel["unit_003"] = createCSTR(1.0);
	jsonModel["unit_004"] = createGeneralRateModel();
	jsonModel["unit_005"] = createOutlet();
	jsonModel["unit_006"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, -1.0, -1.0, 1.0,
	                                             1.0, 3.0, -1.0, -1.0, 3.0,
	                                             3.0, 4.0, -1.0, -1.0, 4.0,
	                                             4.0, 5.0, -1.0, -1.0, 2.0,
	                                             4.0, 6.0, -1.0, -1.0, 2.0
	                                           });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());

	const mpfr::mpreal abscissa(2);
	const mpfr::mpreal T(3);
	const mpfr::mpreal nTerms(97);

	mpfr::mpreal ref(0);
	{
		const mpfr::mpreal M = (constantFromUnit(jsonModel["unit_000"], T, 0, abscissa, 1) * constantFromUnit(jsonModel["unit_002"], T, 1, abscissa, 1) * 0.25 + constantFromUnit(jsonModel["unit_001"], T, 0, abscissa, 1) * 0.75) * constantFromUnit(jsonModel["unit_003"], T, 4, abscissa, 1);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_004"], T, 4, abscissa, M, nTerms));
	}
	{
		const mpfr::mpreal M = (constantFromUnit(jsonModel["unit_000"], T, 0, abscissa, 1) * constantFromUnit(jsonModel["unit_002"], T, 1, abscissa, 1) * 0.25 + constantFromUnit(jsonModel["unit_001"], T, 0, abscissa, 1) * 0.75);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_003"], T, 4, abscissa, M, nTerms));
	}
	{
		const mpfr::mpreal M = constantFromUnit(jsonModel["unit_000"], T, 0, abscissa, 1);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_002"], T, 1, abscissa, M, nTerms));
	}

	const mpfr::mpreal truncErr = sys->truncationError(abscissa, T, nTerms);
	CHECK(abs((truncErr - ref) / ref) <= mpfr::mpreal("1e-16"));
}

TEST_CASE("Truncation error DAG GRM-CSTR", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(20));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 7;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);
	jsonModel["unit_002"] = createGeneralRateModel();
	jsonModel["unit_003"] = createCSTR(1.0);
	jsonModel["unit_004"] = createGeneralRateModel();
	jsonModel["unit_005"] = createOutlet();
	jsonModel["unit_006"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, -1.0, -1.0, 1.0,
	                                             1.0, 3.0, -1.0, -1.0, 3.0,
	                                             3.0, 4.0, -1.0, -1.0, 4.0,
	                                             4.0, 5.0, -1.0, -1.0, 2.0,
	                                             4.0, 6.0, -1.0, -1.0, 2.0
	                                           });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());

	const mpfr::mpreal abscissa(2);
	const mpfr::mpreal T(3);
	const mpfr::mpreal nTerms(97);

	mpfr::mpreal ref(0);
	{
		const mpfr::mpreal M = (constantFromUnit(jsonModel["unit_000"], T, 0, abscissa, 1) * constantFromUnit(jsonModel["unit_002"], T, 1, abscissa, 1) * 0.25 + constantFromUnit(jsonModel["unit_001"], T, 0, abscissa, 1) * 0.75) * constantFromUnit(jsonModel["unit_003"], T, 4, abscissa, 1);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_004"], T, 4, abscissa, M, nTerms));
	}
	{
		const mpfr::mpreal M = (constantFromUnit(jsonModel["unit_000"], T, 0, abscissa, 1) * constantFromUnit(jsonModel["unit_002"], T, 1, abscissa, 1) * 0.25 + constantFromUnit(jsonModel["unit_001"], T, 0, abscissa, 1) * 0.75);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_003"], T, 4, abscissa, M, nTerms));
	}
	{
		const mpfr::mpreal M = constantFromUnit(jsonModel["unit_000"], T, 0, abscissa, 1);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_002"], T, 1, abscissa, M, nTerms));
	}

	const mpfr::mpreal truncErr = sys->truncationError(abscissa, T, nTerms);
	CHECK(abs((truncErr - ref) / ref) <= mpfr::mpreal("1e-16"));
}

TEST_CASE("Truncation error two chains", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(20));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 7;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);
	jsonModel["unit_002"] = createGeneralRateModel();
	jsonModel["unit_003"] = createCSTR(1.0);
	jsonModel["unit_004"] = createGeneralRateModel();
	jsonModel["unit_005"] = createOutlet();
	jsonModel["unit_006"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, -1.0, -1.0, 1.0,
	                                             3.0, 5.0, -1.0, -1.0, 1.0,
	                                             1.0, 4.0, -1.0, -1.0, 2.0,
	                                             4.0, 6.0, -1.0, -1.0, 2.0
	                                           });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());

	const mpfr::mpreal abscissa(2);
	const mpfr::mpreal T(3);
	const mpfr::mpreal nTerms(97);

	mpfr::mpreal ref(0);
	{
		const mpfr::mpreal M = constantFromUnit(jsonModel["unit_000"], T, 0, abscissa, 1) * constantFromUnit(jsonModel["unit_002"], T, 1, abscissa, 1);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_003"], T, 1, abscissa, M, nTerms));
	}
	{
		const mpfr::mpreal M = constantFromUnit(jsonModel["unit_000"], T, 0, abscissa, 1);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_002"], T, 1, abscissa, M, nTerms));
	}
	{
		const mpfr::mpreal M = constantFromUnit(jsonModel["unit_001"], T, 0, abscissa, 1);
		ref = max(ref, truncationErrorFromUnit(jsonModel["unit_004"], T, 2, abscissa, M, nTerms));
	}

	const mpfr::mpreal truncErr = sys->truncationError(abscissa, T, nTerms);
	CHECK(abs((truncErr - ref) / ref) <= mpfr::mpreal("1e-16"));
}

TEST_CASE("Error estimate DAG GRM", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 7;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);
	jsonModel["unit_002"] = createGeneralRateModel(0.1);
	jsonModel["unit_003"] = createGeneralRateModel(1.0);
	jsonModel["unit_004"] = createGeneralRateModel();
	jsonModel["unit_005"] = createOutlet();
	jsonModel["unit_006"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, -1.0, -1.0, 1.0,
	                                             1.0, 3.0, -1.0, -1.0, 3.0,
	                                             3.0, 4.0, -1.0, -1.0, 4.0,
	                                             4.0, 5.0, -1.0, -1.0, 2.0,
	                                             4.0, 6.0, -1.0, -1.0, 2.0
	                                           });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());

	const mpfr::mpreal error("1e-20");
	const mpfr::mpreal weight("0.5");
	const mpfr::mpreal T = secTimes.back() / 2 * mpfr::mpreal(1.01);
	const mpfr::mpreal abscissaSafety = 0.25;
	const mpfr::mpreal abscissa = log1p(sys->timeDomainUpperBound() / (weight * error)) / (2 * T) + abscissaSafety;
	const mpfr::mpreal consError = sys->timeDomainUpperBound() / expm1(2 * T * abscissa);

	const mpfr::mpreal nTerms = ceil(sys->inverseTruncationError(abscissa, T, (1 - weight) * error * T / exp(2 * abscissa * T)));
	const mpfr::mpreal truncErr = sys->truncationError(abscissa, T, nTerms);
	CHECK(consError + exp(2 * abscissa * T) / T * truncErr <= error);
}

TEST_CASE("Error estimate DAG CSTR", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 7;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);
	jsonModel["unit_002"] = createCSTR(0.1);
	jsonModel["unit_003"] = createCSTR(1.0);
	jsonModel["unit_004"] = createGeneralRateModel();
	jsonModel["unit_005"] = createOutlet();
	jsonModel["unit_006"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, -1.0, -1.0, 1.0,
	                                             1.0, 3.0, -1.0, -1.0, 3.0,
	                                             3.0, 4.0, -1.0, -1.0, 4.0,
	                                             4.0, 5.0, -1.0, -1.0, 2.0,
	                                             4.0, 6.0, -1.0, -1.0, 2.0
	                                           });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());

	const mpfr::mpreal error("1e-20");
	const mpfr::mpreal weight("0.5");
	const mpfr::mpreal T = secTimes.back() / 2 * mpfr::mpreal(1.01);
	const mpfr::mpreal abscissaSafety = 0.25;
	const mpfr::mpreal abscissa = log1p(sys->timeDomainUpperBound() / (weight * error)) / (2 * T) + abscissaSafety;
	const mpfr::mpreal consError = sys->timeDomainUpperBound() / expm1(2 * T * abscissa);

	const mpfr::mpreal nTerms = ceil(sys->inverseTruncationError(abscissa, T, (1 - weight) * error * T / exp(2 * abscissa * T)));
	const mpfr::mpreal truncErr = sys->truncationError(abscissa, T, nTerms);
	CHECK(consError + exp(2 * abscissa * T) / T * truncErr <= error);
}

TEST_CASE("Error estimate DAG GRM-CSTR", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 7;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);
	jsonModel["unit_002"] = createGeneralRateModel();
	jsonModel["unit_003"] = createCSTR(1.0);
	jsonModel["unit_004"] = createGeneralRateModel();
	jsonModel["unit_005"] = createOutlet();
	jsonModel["unit_006"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, -1.0, -1.0, 1.0,
	                                             1.0, 3.0, -1.0, -1.0, 3.0,
	                                             3.0, 4.0, -1.0, -1.0, 4.0,
	                                             4.0, 5.0, -1.0, -1.0, 2.0,
	                                             4.0, 6.0, -1.0, -1.0, 2.0
	                                           });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());

	const mpfr::mpreal error("1e-20");
	const mpfr::mpreal weight("0.5");
	const mpfr::mpreal T = secTimes.back() / 2 * mpfr::mpreal(1.01);
	const mpfr::mpreal abscissaSafety = 0.25;
	const mpfr::mpreal abscissa = log1p(sys->timeDomainUpperBound() / (weight * error)) / (2 * T) + abscissaSafety;
	const mpfr::mpreal consError = sys->timeDomainUpperBound() / expm1(2 * T * abscissa);

	const mpfr::mpreal nTerms = ceil(sys->inverseTruncationError(abscissa, T, (1 - weight) * error * T / exp(2 * abscissa * T)));
	const mpfr::mpreal truncErr = sys->truncationError(abscissa, T, nTerms);
	CHECK(consError + exp(2 * abscissa * T) / T * truncErr <= error);
}

TEST_CASE("Error estimate two chains", "[Estimate]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 7;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);
	jsonModel["unit_002"] = createGeneralRateModel();
	jsonModel["unit_003"] = createCSTR(1.0);
	jsonModel["unit_004"] = createGeneralRateModel();
	jsonModel["unit_005"] = createOutlet();
	jsonModel["unit_006"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, -1.0, -1.0, 1.0,
	                                             3.0, 5.0, -1.0, -1.0, 1.0,
	                                             1.0, 4.0, -1.0, -1.0, 2.0,
	                                             4.0, 6.0, -1.0, -1.0, 2.0
	                                           });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::io::JsonParameterProvider jpp(jsonModel);

	casema::ModelBuilder mb;
	casema::model::ModelSystem* const sys = mb.createSystem(jpp);
	sys->setSectionTimes(secTimes.data(), secTimes.size());

	REQUIRE(sys->hasValidEstimate());

	const mpfr::mpreal error("1e-20");
	const mpfr::mpreal weight("0.5");
	const mpfr::mpreal T = secTimes.back() / 2 * mpfr::mpreal(1.01);
	const mpfr::mpreal abscissaSafety = 0.25;
	const mpfr::mpreal abscissa = log1p(sys->timeDomainUpperBound() / (weight * error)) / (2 * T) + abscissaSafety;
	const mpfr::mpreal consError = sys->timeDomainUpperBound() / expm1(2 * T * abscissa);

	const mpfr::mpreal nTerms = ceil(sys->inverseTruncationError(abscissa, T, (1 - weight) * error * T / exp(2 * abscissa * T)));
	const mpfr::mpreal truncErr = sys->truncationError(abscissa, T, nTerms);
	CHECK(consError + exp(2 * abscissa * T) / T * truncErr <= error);
}
