// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>
#include <boost/uuid/uuid_io.hpp>
#include <catch2/catch.hpp>

#include "SimulationMaster.h"
#include "lb/BuildSystemInterface.h"
#include "lb/lattices/D3Q19.h"
#include "Traits.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "tests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    TEST_CASE_METHOD(helpers::FolderTestFixture,
		     "SadCellIntegrationTests", "[redblood][.long]") {
      typedef Traits<>::Reinstantiate<lb::lattices::D3Q19, lb::GuoForcingLBGK>::Type Traits;
      typedef hemelb::redblood::CellController<Traits> CellControl;
      typedef SimulationMaster<Traits> MasterSim;

      redblood::KruegerMeshIO msh_io = {};
      redblood::VTKMeshIO vtk_io = {};

      CopyResourceToTempdir("large_cylinder_rbc.xml");
      CopyResourceToTempdir("large_cylinder.gmy");
      CopyResourceToTempdir("rbc_ico_1280.msh");

      ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 10000);
      ModifyXMLInput("large_cylinder_rbc.xml", { "redbloodcells",
	    "cells",
	    "cell",
	    "shape",
	    "mesh_path" },
	"rbc_ico_1280.msh");
      ModifyXMLInput("large_cylinder_rbc.xml", { "inlets",
	    "inlet",
	    "condition",
	    "mean",
	    "value" },
	0);
      DeleteXMLInput("large_cylinder_rbc.xml", { "redbloodcells", "cells" });
      DeleteXMLInput("large_cylinder_rbc.xml", { "inlets", "inlet", "flowextension" });
      DeleteXMLInput("large_cylinder_rbc.xml", { "outlets", "outlet", "flowextension" });

      constexpr int argc = 7;
      char const * argv[argc] = {
	"hemelb",
	"-in",
	"large_cylinder_rbc.xml",
	"-i",
	"0",
	"-ss",
	"1111"
      };
      hemelb::configuration::CommandLine options(argc, argv);

      auto master = std::make_shared<MasterSim>(options, Comms());

      SECTION("integration test") {
	auto const normal = msh_io.readFile(resources::Resource("rbc_ico_1280.msh").Path(), true);
	auto const cell = std::make_shared<redblood::Cell>(normal->vertices, normal);
	auto const & converter = master->GetUnitConverter();
	auto const deformed = msh_io.readFile(resources::Resource("sad.msh").Path(), true);
	auto const sadcell = std::make_shared<redblood::Cell>(deformed->vertices, normal);
	auto const scale = converter.ConvertToLatticeUnits("m", 4e-6);
	// scale and positions cell in the middle somewhere
	cell->SetScale(scale);
	*cell *= scale;
	sadcell->SetScale(scale);
	*sadcell *= 1e0 / converter.GetVoxelSize();
	*sadcell += converter.ConvertPositionToLatticeUnits(PhysicalPosition(0, 0, 0))
	  - sadcell->GetBarycenter();
	*cell += converter.ConvertPositionToLatticeUnits(PhysicalPosition(0, 0, 0))
	  - cell->GetBarycenter();
	vtk_io.writeFile("/tmp/ideal.vtp", *cell, converter);
	vtk_io.writeFile("/tmp/deformed.vtp", *sadcell, converter);
	sadcell->moduli.bending = 0.0000375;
	sadcell->moduli.surface = 1e0;
	sadcell->moduli.volume = 1e0;
	sadcell->moduli.dilation = 0.5;
	sadcell->moduli.strain = 0.0006;

	auto controller = std::static_pointer_cast<CellControl>(master->GetCellController());
	controller->AddCell(sadcell);
	std::vector<PhysicalEnergy> energies;
	energies.push_back( (*sadcell)());
	controller->AddCellChangeListener([&energies](CellContainer const &cells) {
	    energies.push_back((**cells.begin())());
	  });

	controller->AddCellChangeListener([&vtk_io, &converter](const hemelb::redblood::CellContainer &cells) {
	    static int iter = 0;
	    for (auto cell : cells) {
	      if(iter % 1000 == 0) {
		std::stringstream filename;
		filename << cell->GetTag() << "_t_" << iter << ".vtp";
		vtk_io.writeFile(filename.str(), *cell, converter);
	      }
	    }
	    ++iter;
	  });

	// run the simulation
	master->RunSimulation();
	vtk_io.writeFile("/tmp/reformed.vtp", *sadcell, converter);

	*sadcell += converter.ConvertPositionToLatticeUnits(PhysicalPosition(0, 0, 0))
	  - sadcell->GetBarycenter();
	vtk_io.writeFile("/tmp/reformed_centered.vtp", *sadcell, converter);

	AssertPresent("results/report.txt");
	AssertPresent("results/report.xml");
      }
    }

  } // namespace tests
} // namespace hemelb