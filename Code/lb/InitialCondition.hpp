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

  const auto num_points = centrelineCoordinate.size();
  
  for (site_t i = 0; i < latDat->GetDomain().GetLocalFluidSiteCount(); i++) {

    // obtain lattice coordinates
    auto site = latDat->GetSite(i);
    const auto coordinates = site.GetGlobalSiteCoords();

    // obtain the closest centreline points
    double closest_distance = std::numeric_limits<double>::max();
    site_t closest_point_index = 0;

    for (site_t j = 1; j < num_points; j++) {
      auto centreline_distance = (centrelineCoordinate[j] - coordinates).GetMagnitude(); 
      if (centreline_distance < closest_distance) {
        closest_distance = centreline_distance;
        closest_point_index = j;
      }
    }

    // use the first and closest points to define the centreline
    auto centreline_point1 = centrelineCoordinate[0];
    auto centreline_point2 = centrelineCoordinate[closest_point_index];

    LatticeDistance radius = radiusDistance[closest_point_index];
    LatticeDistance length = (centreline_point1 - centreline_point2).GetMagnitude(); 

    auto velocity = velocityCoordinate[closest_point_index]; 

    auto pressure_point1 = pressureMagnitude[0];
    auto pressure_point2 = pressureMagnitude[closest_point_index]; // issue p_diff = 0

    // extrapolate from centreline properties to lattice properties
    auto unit_vector = (centreline_point2 - centreline_point1).GetNormalised();
    auto centreline_lattice_vector = coordinates - centreline_point1;
    auto projection_distance = Dot(centreline_lattice_vector, unit_vector);
    auto perpendicular_vector = centreline_lattice_vector - projection_distance * unit_vector;
    auto perpendicular_distance = perpendicular_vector.GetMagnitude();

    LatticePressure lattice_pressure = pressure_point1 - (pressure_point1 - pressure_point2) * projection_distance / length;
    LatticeVelocity lattice_velocity = (velocity * unit_vector) * (1 - (perpendicular_distance * perpendicular_distance) / (radius * radius));

    // convert velocity and pressure to density and momentum
    LatticeDensity density = lattice_pressure / Cs2;
    LatticeMomentum lattice_momentum = density * lattice_velocity;

    // std::cout << ioComms.Rank() << "|" << coordinates[0] << "|" << coordinates[1] << "|" << coordinates[2] 
    // << "|" << lattice_pressure << "|" << lattice_velocity[0] << "|" << lattice_velocity[1] << "|" << lattice_velocity[2] 
    // << "|" << perpendicular_distance  << "|" << radius << std::endl;
    
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
