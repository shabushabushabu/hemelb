// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_EXTRACTION_DUMMYDATASOURCE_H
#define HEMELB_TESTS_EXTRACTION_DUMMYDATASOURCE_H

#include <vector>

#include "util/Vector3D.h"
#include "extraction/IterableDataSource.h"
#include "util/Matrix3D.h"

#include "tests/helpers/RandomSource.h"

namespace hemelb::tests
{
    class DummyDataSource : public extraction::IterableDataSource
    {
    public:
          DummyDataSource() :
              randomNumberGenerator(1358), siteCount(64), location(0), gridPositions(siteCount),
                  pressures(siteCount), velocities(siteCount), voxelSize(0.3e-3),
                  origin(0.034, 0.001, 0.074)
          {
            unsigned ijk = 0;

            // Set up 3d index of sites
            for (unsigned i = 0; i < 4; ++i)
            {
              for (unsigned j = 0; j < 4; ++j)
              {
                for (unsigned k = 0; k < 4; ++k)
                {
                  gridPositions[ijk] = {i, j, k};
                  ++ijk;
                }
              }
            }

          }

          void FillFields()
          {
            // Fill in pressure & velocities as in the Python
            for (unsigned i = 0; i < siteCount; ++i)
            {
              pressures[i] = 80. + 2. * randomNumberGenerator.uniform();
              for (unsigned j = 0; j < 3; ++j)
              {
                velocities[i][j] = 0.01 * randomNumberGenerator.uniform();
              }
            }

          }

          void Reset()
          {
            location = 0 - 1;
          }

          bool ReadNext()
          {
            ++location;
            return location < siteCount;
          }

          hemelb::util::Vector3D<site_t> GetPosition() const
          {
            return gridPositions[location];
          }

          distribn_t GetVoxelSize() const override
          {
            return voxelSize;
          }

          const util::Vector3D<distribn_t>& GetOrigin() const override
          {
            return origin;
          }
          hemelb::extraction::FloatingType GetPressure() const override
          {
            return pressures[location];
          }
          hemelb::util::Vector3D<hemelb::extraction::FloatingType> GetVelocity() const override
          {
            return velocities[location];
          }
          hemelb::extraction::FloatingType GetShearStress() const override
          {
            return 0.f;
          }
          hemelb::extraction::FloatingType GetVonMisesStress() const override
          {
            return 0.f;
          }
          hemelb::extraction::FloatingType GetShearRate() const override
          {
            return 0.f;
          }
          util::Matrix3D GetStressTensor() const override
          {
            //! @todo: #177 add constructor with initialisation to Matrix3D
            util::Matrix3D retValue;
            return retValue;
          }

          util::Vector3D<PhysicalStress> GetTraction() const override
          {
            util::Vector3D<PhysicalStress> retValue(0);
            return retValue;
          }

          util::Vector3D<PhysicalStress> GetTangentialProjectionTraction() const override
          {
            util::Vector3D<PhysicalStress> retValue(0);
            return retValue;
          }

          distribn_t const* GetDistribution() const override
          {
            double fval = 1.0;
            double *p_fval = &fval;
            return p_fval;
          }

	  unsigned GetNumVectors() const override
	  {
	    return 15;
	  }

          bool IsValidLatticeSite(const hemelb::util::Vector3D<site_t>&) const override
          {
            return true;
          }
          bool IsAvailable(const hemelb::util::Vector3D<site_t>&) const override
          {
            return true;
          }
          bool IsWallSite(const util::Vector3D<site_t>& location) const override
          {
            /// @todo: #375 This method is not covered by any test. I'm leaving it unimplemented.
            assert(false);

            return false;
          }

    private:
          helpers::RandomSource randomNumberGenerator;
          const site_t siteCount;
          site_t location;
          std::vector<hemelb::util::Vector3D<site_t> > gridPositions;
          std::vector<distribn_t> pressures;
          std::vector<hemelb::util::Vector3D<hemelb::extraction::FloatingType> > velocities;
          distribn_t voxelSize;
          hemelb::util::Vector3D<distribn_t> origin;
    };

}
#endif // HEMELB_TESTS_EXTRACTION_DUMMYDATASOURCE_H
