// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_REDBLOOD_PARALLEL_GRAPHBASEDCOMMUNICATION_H
#define HEMELB_REDBLOOD_PARALLEL_GRAPHBASEDCOMMUNICATION_H

#include "geometry/Domain.h"
#include "geometry/LookupTree.h"
#include "redblood/types.h"
#include "reporting/timers_fwd.h"
namespace hemelb::net { class MpiCommunicator; }

namespace hemelb::redblood::parallel
{

    struct VectorOctreeOrdering {
        constexpr bool operator()( const Vec16& lhs, const Vec16& rhs ) const {
            return geometry::octree::ijk_to_oct(lhs) < geometry::octree::ijk_to_oct(rhs);
        }
    };

    //! @todo #668 Should this be a std::unordered_map instead? The
    //! map is to be created once (and never modified again) and
    //! look-up heavy later on. Possibly a boost::flatmap/sorted
    //! vector might be more performant.
    using GlobalCoordsToProcMap = std::map<Vec16, proc_t, VectorOctreeOrdering>;

      /**
       * Computes a map that allows looking up the process owning certain lattice sites (by global lattice coordinate).
       * The map is restricted to the local lattices sites and those owned by the neighbours defined in the input graph MPI communicator.
       * @param comm graph MPI communicator
       * @param domain object with distributed information about the lattice
       * @return map between lattice coordinates and process owning that lattice (restricted to graph neighbours)
       */
      GlobalCoordsToProcMap ComputeGlobalCoordsToProcMap(net::MpiCommunicator const &comm,
                                                         const geometry::Domain &domain);

      //! \brief Compute neighbourhood based on checking the minimum distance between every pair of subdomains and declaring them neighbours if this is shorter than the RBCs effective size.
      std::vector<int> ComputeProcessorNeighbourhood(net::MpiCommunicator const &comm,
                                                                  geometry::Domain &domain,
                                                                  LatticeDistance cellsEffectiveSize);

      // In order to compute the graph neighbourhood we assume that cells elongate maximum MAXIMUM_SIZE_TO_RADIUS_RATIO times the radius
      static const LatticeDistance MAXIMUM_SIZE_TO_RADIUS_RATIO = 5.0;

      LatticeDistance ComputeCellsEffectiveSize(TemplateCellContainer const& cellTemplates);

      //! \brief Generates a graph communicator describing the data dependencies for interpolation and spreading
      net::MpiCommunicator CreateGraphComm(net::MpiCommunicator const &comm,
                                           geometry::Domain &domain,
                                           const std::shared_ptr<TemplateCellContainer>& cellTemplates,
                                           hemelb::reporting::Timers &timings);
}
#endif
