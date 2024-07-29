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
  // LatticeType::CalculateFeq(1.0, 0.0, 0.0, 0.0, f_eq);

  // read centerline object
  // TODO: change to list of centerline coordinates, convert to lattice unit (builder)
  const auto centerline_point1 = centrelineCoordinate[0];
  const auto centerline_point2 = centrelineCoordinate.back();

  const float pressure_point1 = pressureMagnitude[0]; // 0.000000;
  const float pressure_point2 = pressureMagnitude.back(); // 0.000000;
  std::cout << "pressure_point1: " << pressure_point1 << std::endl; // REMOVE
  std::cout << "pressure_point2: " << pressure_point2 << std::endl; // REMOVE

  const auto velocity_point1 = velocityCoordinate[0];
  const auto velocity_point2 = velocityCoordinate.back();

  const float radius_point1 = radiusDistance[0];
  const float radius_point2 = radiusDistance.back();

  // TODO: Remove hardcode. Change to read from configuration
  auto length = (centerline_point1 - centerline_point2).GetMagnitude();

  // TODO: Move inside the loop after the correct centerline point is obtained
  const auto unit_vector = (centerline_point2 - centerline_point1).GetNormalised();
  // std::cout << "unit_vec" << unit_vector[0] << unit_vector[1] << unit_vector[2] << std::endl; // REMOVE
  
  for (site_t i = 0; i < latDat->GetDomain().GetLocalFluidSiteCount(); i++) {

    // obtain lattice coordinates
    auto site = latDat->GetSite(i);
    const auto coordinates = site.GetGlobalSiteCoords();

    // TODO: Add calculation to obtain the closest centerline points (currently use the further end)

    // extrapolate from centerline properties to lattice properties
    const auto centerline_lattice_vector = coordinates - centerline_point1;
    auto projection_distance = Dot(centerline_lattice_vector, unit_vector);
    const auto  perpendicular_vector = centerline_lattice_vector - projection_distance * unit_vector;
    auto perpendicular_distance = perpendicular_vector.GetMagnitude();

    LatticePressure lattice_pressure = pressure_point1 - (pressure_point1 - pressure_point2) * projection_distance / length;
    LatticeVelocity lattice_velocity = velocity_point1 * (1 - (perpendicular_distance * perpendicular_distance) / (radius_point1 * radius_point1));

    // convert velocity and pressure to density and momentum
    LatticeDensity density = lattice_pressure / Cs2;
    LatticeMomentum lattice_momentum; // TODO: Simplify
    // density * lattice_velocity
    lattice_momentum = density * lattice_velocity;
    // lattice_momentum.x() = density * lattice_velocity.x();
    // lattice_momentum.y() = density * lattice_velocity.y();
    // lattice_momentum.z() = density * lattice_velocity.z();
    
    LatticeType::CalculateFeq(density, lattice_momentum, f_eq);
    // LatticeType::CalculateFeq(density, lattice_momentum.x(), lattice_momentum.y(), lattice_momentum.z() , f_eq);

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
