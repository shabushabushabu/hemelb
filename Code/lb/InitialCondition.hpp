// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_INITIALCONDITION_HPP
#define HEMELB_LB_INITIALCONDITION_HPP

#include "lb/InitialCondition.h"
#include "extraction/LocalDistributionInput.h"

namespace hemelb {
  namespace lb {
    template<class LatticeType>
    struct FSetter {
      using result_type = void;
      void operator()(std::monostate) {
	throw Exception() << "No initial condition specified";
      }

      template <typename T>
      void operator()(T&& t) const {
	t.template SetFs<LatticeType>(latDat, ioComms);
      }
      geometry::FieldData* latDat;
      const net::IOCommunicator& ioComms;
    };

    
    template<class LatticeType>
    void InitialCondition::SetFs(geometry::FieldData* latDat, const net::IOCommunicator& ioComms) const {
      // Since as far as std::visit is concerned, `this` is not a
      // variant.
      auto this_var = static_cast<ICVar const*>(this);

      std::visit(FSetter<LatticeType>{latDat, ioComms}, *this_var);
    }

template<class LatticeType>
void EquilibriumInitialCondition::SetFs(geometry::FieldData* latDat, const net::IOCommunicator& ioComms) const {
  distribn_t f_eq[LatticeType::NUMVECTORS];
  LatticeType::CalculateFeq(density, mom_x, mom_y, mom_z, f_eq);

  for (site_t i = 0; i < latDat->GetDomain().GetLocalFluidSiteCount(); i++) {
    distribn_t* f_old_p = this->GetFOld(latDat, i * LatticeType::NUMVECTORS);
    distribn_t* f_new_p = this->GetFNew(latDat, i * LatticeType::NUMVECTORS);

    for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++) {
      f_new_p[l] = f_old_p[l] = f_eq[l];
    }
  }

}

template<class LatticeType>
void CentrelineInitialCondition::SetFs(geometry::FieldData* latDat, const net::IOCommunicator& ioComms) const {
  distribn_t f_eq[LatticeType::NUMVECTORS];
  LatticeType::CalculateFeq(density, mom_x, mom_y, mom_z, f_eq);

  // read centerline object
  // TODO: change to list of centerline coordinates, convert to lattice unit (builder)
  const util::Vector3D<float> centerline_point1 = {0.000000,-0.000000,50.061935};
  const util::Vector3D<float> centerline_point2 = {0.000000,0.000000,-37.801815};

  const float pressure_point1 = 2.666440; // scalar
  const float pressure_point2 = 0.323606;

  const util::Vector3D<float> velocity_point1 = {0.000000,0.000000,0.266644};
  const util::Vector3D<float> velocity_point2 = {0.000000,0.000000,0.266644};

  // TODO: Remove hardcode. Change to read from configuration
  auto radius = 0.5;
  auto length = (centerline_point1 - centerline_point2).GetMagnitude();
  configuration::FluidInfo fluidInfo;
  double fluid_viscosity = fluidInfo.viscosity_Pas; // need to convert to lattice unit?

  // TODO: Move inside the loop after the correct centerline point is obtain
  const util::Vector3D<float> unit_vector = (centerline_point2 - centerline_point1).GetNormalised();
  // std::cout << "unit_vec" << unit_vector[0] << unit_vector[1] << unit_vector[2] << std::endl; // REMOVE
  
  for (site_t i = 0; i < latDat->GetDomain().GetLocalFluidSiteCount(); i++) {

    // obtain lattice coordinates
    auto site = latDat->GetSite(i);
    const util::Vector3D<site_t> coordinates = site.GetGlobalSiteCoords();

    // TODO: Add calculation to obtain the closest centerline points (currently use the further end)

    // extrapolate from centerline properties to lattice properties
    const util::Vector3D<float> centerline_lattice_vector = coordinates - centerline_point1;
    auto projection_distance = Dot(centerline_lattice_vector, unit_vector);
    const util::Vector3D<float>  perpendicular_vector = centerline_lattice_vector - projection_distance * unit_vector;
    auto perpendicular_distance = perpendicular_vector.GetMagnitude(); 
    // fix scalar P
    auto lattice_pressure = pressure_point1 - (pressure_point1 - pressure_point2) * projection_distance / length;
    // LatticeVelocity lattice_velocity = (pressure_point1 - pressure_point2) / (4 * fluid_viscosity * length) * (radius * radius - perpendicular_distance * perpendicular_distance);
    LatticeVelocity lattice_velocity = velocity_point1 * (1 - (perpendicular_distance * perpendicular_distance) / (radius * radius));

    // conversions TBC
    LatticeDensity density = lattice_pressure / Cs2;
    // distribn_t mom_x = static_cast<distribn_t>(density * lattice_velocity[0]);
    // distribn_t mom_y = static_cast<distribn_t>(density * lattice_velocity[1]);
    // distribn_t mom_z = static_cast<distribn_t>(density * lattice_velocity[2]);

    LatticeMomentum lattice_momentum;
    lattice_momentum.x() = density * lattice_velocity.x();
    lattice_momentum.y() = density * lattice_velocity.y();
    lattice_momentum.z() = density * lattice_velocity.z();

    // LatticeMomentum lattice_momentum = {density * lattice_velocity[0], density * lattice_velocity[1], density * lattice_velocity[2]};

    std::cout << "mom_x: " << lattice_momentum.x() << " mom_y: " << lattice_momentum.y() << " mom_z: " << lattice_momentum.z() << std::endl;

    
    // should accept LatticeMomentum 3D
    LatticeType::CalculateFeq(density, lattice_momentum, f_eq);

    // get dist 
    distribn_t* f_old_p = this->GetFOld(latDat, i * LatticeType::NUMVECTORS);
    distribn_t* f_new_p = this->GetFNew(latDat, i * LatticeType::NUMVECTORS);

    for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++) {
      f_new_p[l] = f_old_p[l] = f_eq[l];
    }
  }
}


    template<class LatticeType>
    void CheckpointInitialCondition::SetFs(geometry::FieldData* latDat, const net::IOCommunicator& ioComms) const {
      auto distributionInputPtr = std::make_unique<extraction::LocalDistributionInput>(cpFile, maybeOffFile, ioComms);
      distributionInputPtr->LoadDistribution(latDat, initial_time);
    }

  }
}

#endif
