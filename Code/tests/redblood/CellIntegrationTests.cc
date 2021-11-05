// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cstdio>

#include <catch2/catch.hpp>

#include "Traits.h"
#include "SimulationMaster.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;
    class CellIntegrationTests : public helpers::FolderTestFixture
    {
      using Traits = Traits<>::ChangeKernel<lb::GuoForcingLBGK>::Type;
      using CellControll = CellController<Traits>;
      using MasterSim = SimulationMaster<Traits>;

    public:
      CellIntegrationTests() : FolderTestFixture(), timings(Comms()) {
	argv[0] = "hemelb";
	argv[1] = "-in";
	argv[2] = "large_cylinder.xml";
	argv[3] = "-i";
	argv[4] = "1";
	argv[5] = "-ss";
	argv[6] = "1111";

	CopyResourceToTempdir("red_blood_cell.txt");
	TiXmlDocument doc(resources::Resource("large_cylinder.xml").Path());
	CopyResourceToTempdir("large_cylinder.xml");
	std::vector<std::string> intel;
	intel.push_back("simulation");
	intel.push_back("steps");
	intel.push_back("value");
	ModifyXMLInput("large_cylinder.xml", std::move(intel), 8);
	CopyResourceToTempdir("large_cylinder.gmy");
	options = std::make_shared<hemelb::configuration::CommandLine>(argc, argv);

	auto cell = std::make_shared<Cell>(icoSphere(4));
	templates = std::make_shared<TemplateCellContainer>();
	(*templates)["icosphere"] = cell;
	cell->moduli = Cell::Moduli(1e-6, 1e-6, 1e-6, 1e-6);
	cells.insert(cell);

	//timings = std::make_unique<reporting::Timers>(Comms());
	master = std::make_shared<MasterSim>(*options, Comms());
	helpers::LatticeDataAccess(&master->GetLatticeData()).ZeroOutForces();
      }

      ~CellIntegrationTests() {
	master->Finalise();
      }

      // No errors when interpolation/spreading hits nodes out of bounds
      void testCellOutOfBounds()
      {
	(*cells.begin())->operator+=(master->GetLatticeData().GetGlobalSiteMins() * 2.0);
	auto controller = std::make_shared<CellControll>(
							 master->GetLatticeData(),
							 cells,
							 templates,
							 timings
							 );
	master->RegisterActor(*controller, 1);
	master->RunSimulation();
	AssertPresent("results/report.txt");
	AssertPresent("results/report.xml");
      }

      // Check that the particles move and result in some force acting on the fluid
      void testIntegration()
      {
	// setup cell position
	auto const &latticeData = master->GetLatticeData();
	auto const mid = LatticePosition(latticeData.GetGlobalSiteMaxes()
					 + latticeData.GetGlobalSiteMins()) * 0.5;
	(**cells.begin()) += mid - (*cells.begin())->GetBarycenter();
	(**cells.begin()) += LatticePosition(0, 0, 8 - mid.z);
	(**cells.begin()) *= 5.0;
	auto controller = std::make_shared<CellControll>(master->GetLatticeData(), cells, templates, timings);
	auto const barycenter = (*cells.begin())->GetBarycenter();

	// run
	master->RegisterActor(*controller, 1);
	master->RunSimulation();

	// check position of cell has changed
	auto const moved = (*cells.begin())->GetBarycenter();
	REQUIRE(Approx(barycenter.x).margin(1e-12) == moved.x);
	REQUIRE(Approx(barycenter.y).margin(1e-12) == moved.y);
	REQUIRE(std::abs(barycenter.z - moved.z) > 1e-8);

	// Check there is force on one of the lattice site near a
	// node node position is guessed at from geometry. This
	// truncates.
	auto const nodepos = LatticeVector{mid + LatticePosition(0, 0, 8 - 5 - mid.z)};
	auto const force = latticeData.GetSite(nodepos).GetForce();
	REQUIRE(std::abs(force.z) > 1e-4);

	AssertPresent("results/report.txt");
	AssertPresent("results/report.xml");
      }

      // Check that the particles move and result in some force acting on the fluid
      void testIntegrationWithoutCells()
      {
	// setup cell position
	CellContainer empty;
	auto empty_tmpl = std::make_shared<TemplateCellContainer>();
	auto controller = std::make_shared<CellControll>(
							 master->GetLatticeData(),
							 empty, empty_tmpl, timings);

	// run
	master->RegisterActor(*controller, 1);
	master->RunSimulation();

	AssertPresent("results/report.txt");
	AssertPresent("results/report.xml");
      }

    private:
      std::shared_ptr<MasterSim> master;
      std::shared_ptr<configuration::CommandLine> options;
      CellContainer cells;
      std::shared_ptr<TemplateCellContainer> templates;
      reporting::Timers timings;
      int const argc = 7;
      char const * argv[7];
    };

    METHOD_AS_TEST_CASE(CellIntegrationTests::testCellOutOfBounds,
			"No errors when interpolation/spreading hits nodes out of bounds",
			"[redblood][.long]");
    METHOD_AS_TEST_CASE(CellIntegrationTests::testIntegration,
			"Check that the particles move and result in some force acting on the fluid",
			"[redblood][.long]");
    METHOD_AS_TEST_CASE(CellIntegrationTests::testIntegrationWithoutCells,
			"Test things work without any cells",
			"[redblood][.long]");

  }
}