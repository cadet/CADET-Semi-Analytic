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
#include "BesselZeros.hpp"

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

	json createGeneralRateModel2D(double radius = 0.01, double colDispersion = 5.75e-8, double colLength = 0.014, double kA = 35.5)
	{
		json j = json::parse(R"json({
			"UNIT_TYPE": "GENERAL_RATE_MODEL_2D",
			"NCOMP": 1,
			"VELOCITY": 5.75e-4,
			"COL_DISPERSION": 5.75e-8,
			"COL_DISPERSION_RADIAL": 1e-6,
			"FILM_DIFFUSION": [6.9e-6],
			"PAR_DIFFUSION": [6.07e-11],
			"PAR_SURFDIFFUSION": [0.0],
			"COL_LENGTH": 0.014,
			"COL_RADIUS": 0.01,
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
				"NRAD": 1,
				"RADIAL_COMPARTMENTS": [],
				"RADIAL_DISC_TYPE": "EQUIDISTANT",
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
		j["COL_RADIUS"] = radius;
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

	void forceJointSolution(json& sys)
	{
		json solver;
		solver["LINEAR_SOLUTION_MODE"] = 1;
		sys["solver"] = solver;
	}
}

TEST_CASE("Joint system solution single port units", "[System]")
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

	casema::ModelBuilder mb;

	casema::io::JsonParameterProvider jppL(jsonModel);
	casema::model::ModelSystem* const sysL = mb.createSystem(jppL);
	sysL->setSectionTimes(secTimes.data(), secTimes.size());

	forceJointSolution(jsonModel);
	casema::io::JsonParameterProvider jppJ(jsonModel);
	casema::model::ModelSystem* const sysJ = mb.createSystem(jppJ);
	sysJ->setSectionTimes(secTimes.data(), secTimes.size());

	std::vector<mpfr::mpcomplex> solL(7);
	std::vector<mpfr::mpcomplex> solJ(7);

	const mpfr::mpcomplex z(1,1);
	casema::model::ModelSystem::Workspace ws = sysL->makeWorkspace();
	sysL->evaluate(z, solL.data(), ws);
	sysJ->evaluate(z, solJ.data(), ws);

	for (int i = 0; i < solL.size(); ++i)
		CHECK(abs(solL[i] - solJ[i]) <= 10 * std::numeric_limits<mpfr::mpreal>::epsilon());
}

TEST_CASE("Joint system solution two chains single port", "[System]")
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

	casema::ModelBuilder mb;

	casema::io::JsonParameterProvider jppL(jsonModel);
	casema::model::ModelSystem* const sysL = mb.createSystem(jppL);
	sysL->setSectionTimes(secTimes.data(), secTimes.size());

	forceJointSolution(jsonModel);
	casema::io::JsonParameterProvider jppJ(jsonModel);
	casema::model::ModelSystem* const sysJ = mb.createSystem(jppJ);
	sysJ->setSectionTimes(secTimes.data(), secTimes.size());

	std::vector<mpfr::mpcomplex> solL(sysL->numOutlets());
	std::vector<mpfr::mpcomplex> solJ(sysL->numOutlets());

	const mpfr::mpcomplex z(1,1);
	casema::model::ModelSystem::Workspace ws = sysL->makeWorkspace();
	sysL->evaluate(z, solL.data(), ws);
	sysJ->evaluate(z, solJ.data(), ws);

	for (int i = 0; i < solL.size(); ++i)
		CHECK(abs(solL[i] - solJ[i]) <= 10 * std::numeric_limits<mpfr::mpreal>::epsilon());
}

TEST_CASE("Joint system solution multiple ports", "[System]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 5;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);

	json grm2D = createGeneralRateModel2D();
	grm2D["COL_RADIUS"] = static_cast<double>(sqrt(1.0 / (grm2D["COL_POROSITY"].get<double>() * grm2D["VELOCITY"].get<double>() * mpfr::const_pi())));
	grm2D["discretization"]["NRAD"] = 2;

	jsonModel["unit_002"] = grm2D;
	jsonModel["unit_003"] = createOutlet();
	jsonModel["unit_004"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, 0.0, 0.0, -1.0, -1.0, 1.0,
	                                             1.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, 0.0, 0.0, -1.0, -1.0, 1.0,
	                                             2.0, 4.0, 1.0, 0.0, -1.0, -1.0, 1.0
	                                           });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::ModelBuilder mb;

	casema::io::JsonParameterProvider jppL(jsonModel);
	casema::model::ModelSystem* const sysL = mb.createSystem(jppL);
	sysL->setSectionTimes(secTimes.data(), secTimes.size());

	forceJointSolution(jsonModel);
	casema::io::JsonParameterProvider jppJ(jsonModel);
	casema::model::ModelSystem* const sysJ = mb.createSystem(jppJ);
	sysJ->setSectionTimes(secTimes.data(), secTimes.size());

	std::vector<mpfr::mpreal> zeros(101);
	casema::besselZerosJ1(zeros.size(), zeros.data());
	sysL->setBesselZeros(zeros.data(), zeros.size());
	sysJ->setBesselZeros(zeros.data(), zeros.size());

	std::vector<mpfr::mpcomplex> solL(sysL->numOutlets());
	std::vector<mpfr::mpcomplex> solJ(sysL->numOutlets());

	const mpfr::mpcomplex z(1,1);
	casema::model::ModelSystem::Workspace ws = sysL->makeWorkspace();
	sysL->evaluate(z, solL.data(), ws);
	sysJ->evaluate(z, solJ.data(), ws);

	for (int i = 0; i < solL.size(); ++i)
		CHECK(abs(solL[i] - solJ[i]) <= 10 * std::numeric_limits<mpfr::mpreal>::epsilon());
}

TEST_CASE("Joint system solution multiple ports single out", "[System]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));

	const double cIn = 2.0;
	const double cIn2 = 3.0;

	json jsonModel;
	jsonModel["NUNITS"] = 4;
	jsonModel["unit_000"] = createInlet(cIn);
	jsonModel["unit_001"] = createInlet(cIn2);

	json grm2D = createGeneralRateModel2D();
	grm2D["COL_RADIUS"] = static_cast<double>(sqrt(1.0 / (grm2D["COL_POROSITY"].get<double>() * grm2D["VELOCITY"].get<double>() * mpfr::const_pi())));
	grm2D["discretization"]["NRAD"] = 2;

	jsonModel["unit_002"] = grm2D;
	jsonModel["unit_003"] = createOutlet();
	jsonModel["connections"] = makeConnections({ 0.0, 2.0, 0.0, 0.0, -1.0, -1.0, 3.0,
	                                             1.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0,
	                                             2.0, 3.0, 0.0, 0.0, -1.0, -1.0, 3.0,
	                                             2.0, 3.0, 1.0, 0.0, -1.0, -1.0, 1.0
	                                           });
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	casema::ModelBuilder mb;

	casema::io::JsonParameterProvider jppL(jsonModel);
	casema::model::ModelSystem* const sysL = mb.createSystem(jppL);
	sysL->setSectionTimes(secTimes.data(), secTimes.size());

	forceJointSolution(jsonModel);
	casema::io::JsonParameterProvider jppJ(jsonModel);
	casema::model::ModelSystem* const sysJ = mb.createSystem(jppJ);
	sysJ->setSectionTimes(secTimes.data(), secTimes.size());

	std::vector<mpfr::mpreal> zeros(101);
	casema::besselZerosJ1(zeros.size(), zeros.data());
	sysL->setBesselZeros(zeros.data(), zeros.size());
	sysJ->setBesselZeros(zeros.data(), zeros.size());

	std::vector<mpfr::mpcomplex> solL(sysL->numOutlets());
	std::vector<mpfr::mpcomplex> solJ(sysL->numOutlets());

	const mpfr::mpcomplex z(1,1);
	casema::model::ModelSystem::Workspace ws = sysL->makeWorkspace();
	sysL->evaluate(z, solL.data(), ws);
	sysJ->evaluate(z, solJ.data(), ws);

	for (int i = 0; i < solL.size(); ++i)
		CHECK(abs(solL[i] - solJ[i]) <= 10 * std::numeric_limits<mpfr::mpreal>::epsilon());
}
