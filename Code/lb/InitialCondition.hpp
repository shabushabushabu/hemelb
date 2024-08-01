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

// template<class LatticeType>
// void CentrelineInitialCondition::SetFs(geometry::FieldData* latDat, const net::IOCommunicator& ioComms) const {
//   distribn_t f_eq[LatticeType::NUMVECTORS];
//   // LatticeType::CalculateFeq(1.0, 0.0, 0.0, 0.0, f_eq);

//   const auto length_centreline = centrelineCoordinate.size();
//   // auto length = (centrelineCoordinate[0] - centrelineCoordinate.back()).GetMagnitude();  
  
//   for (site_t i = 0; i < latDat->GetDomain().GetLocalFluidSiteCount(); i++) {

//     // obtain lattice coordinates
//     auto site = latDat->GetSite(i);
//     const auto coordinates = site.GetGlobalSiteCoords();

//     // obtain the closest centreline points (1 is closer)
//     double closest_distance1 = std::numeric_limits<double>::max();
//     double closest_distance2 = std::numeric_limits<double>::max();
//     site_t closest_point_index1 = 0;
//     site_t closest_point_index2 = length_centreline - 1;
//     for (site_t j = 1; j < length_centreline; j++) {
//       auto centreline_distance = (centrelineCoordinate[j] - coordinates).GetMagnitude(); 
//       if (centreline_distance < closest_distance1) {
//         closest_distance2 = closest_distance1;
//         closest_point_index2 = closest_point_index1;
//         closest_distance1 = centreline_distance;
//         closest_point_index1 = j;
//       } else if (centreline_distance < closest_distance2) {
//         closest_distance2 = centreline_distance;
//         closest_point_index2 = j;
//       }
//     }
//     auto centreline_point1 = centrelineCoordinate[0];
//     auto centreline_point2 = centrelineCoordinate[closest_point_index1];

//     LatticeDistance radius = radiusDistance[closest_point_index1];
//     LatticeDistance length = (centreline_point1 - centreline_point2).GetMagnitude();  

//     auto pressure_point1 = pressureMagnitude[0];
//     auto pressure_point2 = pressureMagnitude[closest_point_index1]; // issue p_diff = 0

//     auto velocity_point1 = velocityCoordinate[0];

//     // extrapolate from centreline properties to lattice properties
//     auto unit_vector = (centreline_point2 - centreline_point1).GetNormalised();
//     auto centreline_lattice_vector = coordinates - centreline_point1;
//     LatticeDistance projection_distance = Dot(centreline_lattice_vector, unit_vector);
//     auto perpendicular_vector = centreline_lattice_vector - projection_distance * unit_vector;
//     LatticeDistance perpendicular_distance = perpendicular_vector.GetMagnitude();

//     LatticePressure lattice_pressure = pressure_point1 - (pressure_point1 - pressure_point2) * projection_distance / length;
//     LatticeVelocity lattice_velocity = velocity_point1 * (1 - (perpendicular_distance * perpendicular_distance) / (radius * radius));

//     // convert velocity and pressure to density and momentum
//     LatticeDensity density = lattice_pressure / Cs2;
//     LatticeMomentum lattice_momentum = density * lattice_velocity;

//     if (i==0) {
//       std::cout << ioComms.Rank() << " pressure point1: " << pressure_point1 << std::endl;
//       std::cout << ioComms.Rank() << " pressure point2: " << pressure_point2 << std::endl;
//       std::cout << ioComms.Rank() << " closest indexes: " << closest_point_index1 << ", " << closest_point_index2 << std::endl;
//       std::cout << ioComms.Rank() << " closest distance: " << closest_distance1 << ", " << closest_distance2 << std::endl;
//       std::cout << ioComms.Rank() << " lattice_pressure: " << lattice_pressure << std::endl;
//       std::cout << ioComms.Rank() << " lattice_velocity: "  << lattice_velocity[0] << lattice_velocity[1] << lattice_velocity[2] << std::endl;
//       std::cout << ioComms.Rank() << " density: " << density << std::endl;
//       std::cout << ioComms.Rank() << " lattice_momentum: " << lattice_momentum[0] << lattice_momentum[1] << lattice_momentum[2] << std::endl;
//     }
    
//     LatticeType::CalculateFeq(density, lattice_momentum, f_eq);

//     // get dist 
//     distribn_t* f_old_p = this->GetFOld(latDat, i * LatticeType::NUMVECTORS);
//     distribn_t* f_new_p = this->GetFNew(latDat, i * LatticeType::NUMVECTORS);

//     for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++) {
//       f_new_p[l] = f_old_p[l] = f_eq[l];
//     }
//   }
// }


template<class LatticeType>
void CentrelineInitialCondition::SetFs(geometry::FieldData* latDat, const net::IOCommunicator& ioComms) const {
  distribn_t f_eq[LatticeType::NUMVECTORS];
  // LatticeType::CalculateFeq(1.0, 0.0, 0.0, 0.0, f_eq);

  // TODO: change to list of centreline coordinates
  const auto centreline_point1 = centrelineCoordinate.back();
  const auto centreline_point2 = centrelineCoordinate[0];

  const float pressure_point1 = pressureMagnitude.back();
  const float pressure_point2 = pressureMagnitude[0];

  const auto velocity_point1 = velocityCoordinate.back();
  const auto velocity_point2 = velocityCoordinate[0];

  const float radius_point1 = radiusDistance.back();
  const float radius_point2 = radiusDistance[0];

  // TODO: Move inside the loop after the correct centreline point is obtained
  auto length = (centreline_point1 - centreline_point2).GetMagnitude();
  const auto unit_vector = (centreline_point2 - centreline_point1).GetNormalised();
  
  for (site_t i = 0; i < latDat->GetDomain().GetLocalFluidSiteCount(); i++) {

    // obtain lattice coordinates
    auto site = latDat->GetSite(i);
    const auto coordinates = site.GetGlobalSiteCoords();

    // TODO: Add calculation to obtain the closest centreline points (currently use the further end)

    // extrapolate from centreline properties to lattice properties
    const auto centreline_lattice_vector = coordinates - centreline_point1;
    auto projection_distance = Dot(centreline_lattice_vector, unit_vector);
    const auto  perpendicular_vector = centreline_lattice_vector - projection_distance * unit_vector;
    auto perpendicular_distance = perpendicular_vector.GetMagnitude();

    LatticePressure lattice_pressure = pressure_point1 - (pressure_point1 - pressure_point2) * projection_distance / length;
    LatticeVelocity lattice_velocity = (velocity_point1 * unit_vector) * (1 - (perpendicular_distance * perpendicular_distance) / (radius_point1 * radius_point1));

    // convert velocity and pressure to density and momentum
    LatticeDensity density = lattice_pressure / Cs2;
    LatticeMomentum lattice_momentum = density * lattice_velocity;

    // std::cout << ioComms.Rank() << "|" << coordinates[0] << "|" << coordinates[1] << "|" << coordinates[2] 
    // << "|" << lattice_pressure << "|" << lattice_velocity[0] << "|" << lattice_velocity[1] << "|" << lattice_velocity[2] 
    // << "|" << perpendicular_distance  << "|" << radius_point1 << std::endl;
    
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
