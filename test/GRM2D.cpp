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

#include <memory>

#include "io/JsonParameterProvider.hpp"
#include "ModelBuilder.hpp"
#include "model/ModelSystem.hpp"
#include "model/UnitOperation.hpp"
#include "model/GeneralRateModel2D.hpp"
#include "model/GeneralRateModel.hpp"
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
			"COL_DISPERSION_RADIAL": 5e-8,
			"FILM_DIFFUSION": [6.9e-6],
			"PAR_DIFFUSION": [6.07e-11],
			"PAR_SURFDIFFUSION": [0.0],
			"COL_LENGTH": 0.014,
			"COL_RADIUS": 0.03,
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
}

TEST_CASE("GRM vs GRM2D single zone", "[Bessel][GRM]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));

	const mpfr::mpreal Q("1e-4");
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	json jsonModel2d = createGeneralRateModel2D();
	jsonModel2d["COL_RADIUS"] = static_cast<double>(sqrt(Q / (jsonModel2d["COL_POROSITY"].get<double>() * jsonModel2d["VELOCITY"].get<double>() * mpfr::const_pi())));

	casema::io::JsonParameterProvider jpp2d(jsonModel2d);

	casema::ModelBuilder mb;
	std::unique_ptr<casema::model::UnitOperation> m2d(mb.createUnitOperation(jpp2d, 0));
	REQUIRE(m2d);
	casema::model::GeneralRateModel2D* const grm2d = reinterpret_cast<casema::model::GeneralRateModel2D*>(m2d.get());

	json jsonModel1d = createGeneralRateModel();
	casema::io::JsonParameterProvider jpp1d(jsonModel1d);
	std::unique_ptr<casema::model::UnitOperation> m1d(mb.createUnitOperation(jpp1d, 1));
	REQUIRE(m1d);
	casema::model::GeneralRateModel* const grm1d = reinterpret_cast<casema::model::GeneralRateModel*>(m1d.get());

	std::vector<mpfr::mpreal> zeros(101);
	casema::besselZerosJ1(zeros.size(), zeros.data());

	m2d->setSectionTimes(secTimes.data(), secTimes.size());
	m2d->setFlowRates(&Q, &Q);
	m2d->setBesselZeros(zeros.data(), zeros.size());

	m1d->setSectionTimes(secTimes.data(), secTimes.size());
	m1d->setFlowRates(&Q, &Q);

	grm2d->velocity(grm1d->velocity());

	Eigen::VectorCmp g1 = Eigen::VectorCmp::Zero(1);
	Eigen::MatrixCmp H1 = Eigen::MatrixCmp::Zero(1, 1);
	Eigen::VectorCmp g2 = Eigen::VectorCmp::Zero(1);
	Eigen::MatrixCmp H2 = Eigen::MatrixCmp::Zero(1, 1);

	const mpfr::mpcomplex z(1, 1);
	m2d->evaluate(z, H2, g2);
	m1d->evaluate(z, H1, g1);

	CHECK((H1-H2).array().abs().maxCoeff() <= 500 * std::numeric_limits<mpfr::mpreal>::epsilon());
	CHECK((g1-g2).array().abs().maxCoeff() <= 500 * std::numeric_limits<mpfr::mpreal>::epsilon());
}

TEST_CASE("GRM vs GRM2D two zones", "[Bessel][GRM]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(30));

	const mpfr::mpreal Q("1e-4");
	const std::vector<mpfr::mpreal> secTimes{0.0, 1.0, 10.0};

	json jsonModel2d = createGeneralRateModel2D();
	jsonModel2d["COL_RADIUS"] = static_cast<double>(sqrt(Q / (jsonModel2d["COL_POROSITY"].get<double>() * jsonModel2d["VELOCITY"].get<double>() * mpfr::const_pi())));
	jsonModel2d["discretization"]["NRAD"] = 3;

	casema::io::JsonParameterProvider jpp2d(jsonModel2d);

	casema::ModelBuilder mb;
	std::unique_ptr<casema::model::UnitOperation> m2d(mb.createUnitOperation(jpp2d, 0));
	REQUIRE(m2d);
	casema::model::GeneralRateModel2D* const grm2d = reinterpret_cast<casema::model::GeneralRateModel2D*>(m2d.get());

	json jsonModel1d = createGeneralRateModel();
	casema::io::JsonParameterProvider jpp1d(jsonModel1d);
	std::unique_ptr<casema::model::UnitOperation> m1d(mb.createUnitOperation(jpp1d, 1));
	REQUIRE(m1d);
	casema::model::GeneralRateModel* const grm1d = reinterpret_cast<casema::model::GeneralRateModel*>(m1d.get());

	std::vector<mpfr::mpreal> zeros(101);
	casema::besselZerosJ1(zeros.size(), zeros.data());

	m2d->setSectionTimes(secTimes.data(), secTimes.size());
	m2d->setBesselZeros(zeros.data(), zeros.size());

	m1d->setSectionTimes(secTimes.data(), secTimes.size());
	m1d->setFlowRates(&Q, &Q);

	grm2d->velocity(grm1d->velocity());

	Eigen::VectorCmp g1 = Eigen::VectorCmp::Zero(1);
	Eigen::MatrixCmp H1 = Eigen::MatrixCmp::Zero(1, 1);
	Eigen::VectorCmp g2 = Eigen::VectorCmp::Zero(3);
	Eigen::MatrixCmp H2 = Eigen::MatrixCmp::Zero(3, 3);

	Eigen::VectorCmp w = Eigen::VectorCmp::Zero(3);
	w(0) = (sqr(grm2d->radialEdges()[1] / grm2d->colRadius()) - sqr(grm2d->radialEdges()[0] / grm2d->colRadius()));
	w(1) = (sqr(grm2d->radialEdges()[2] / grm2d->colRadius()) - sqr(grm2d->radialEdges()[1] / grm2d->colRadius()));
	w(2) = (sqr(grm2d->radialEdges()[3] / grm2d->colRadius()) - sqr(grm2d->radialEdges()[2] / grm2d->colRadius()));

	const mpfr::mpcomplex z(1, 1);
	m2d->evaluate(z, H2, g2);
	m1d->evaluate(z, H1, g1);

	CHECK(abs(H1(0,0)-w.dot(H2.rowwise().sum())) <= 500 * std::numeric_limits<mpfr::mpreal>::epsilon());
	CHECK(abs(g1(0,0)-g2.sum()) <= 500 * std::numeric_limits<mpfr::mpreal>::epsilon());
}
