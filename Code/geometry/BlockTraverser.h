// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_BLOCKTRAVERSER_H
#define HEMELB_GEOMETRY_BLOCKTRAVERSER_H

#include "lb/lattices/D3Q27.h"
#include "geometry/Block.h"
#include "geometry/VolumeTraverser.h"
#include "geometry/SiteTraverser.h"

namespace hemelb::geometry
{
    /**
     *BlockTraverser is used to traverse the blocks in a lattice sequentially.
     */
    class BlockTraverser : public VolumeTraverser
    {
      public:
        BlockTraverser(const Domain& iLatDat);

        /**
         * @override Of the default destructor in VolumeTraverser.
         */
        ~BlockTraverser() override = default;

        [[nodiscard]] site_t CurrentBlockNumber() const;

        util::Vector3D<site_t> GetSiteCoordinatesOfLowestSiteInCurrentBlock();

        const Block& GetCurrentBlockData();

        const Block& GetBlockDataForLocation(const Vec16& iLocation);

        site_t GetBlockSize();

        SiteTraverser GetSiteTraverser();

        /**
         * Gets the number of blocks in the x direction.
         * @override of the abstract function in VolumeTraverser.
         * @return
         */
        [[nodiscard]] U16 GetXCount() const override;

        /**
         * Gets the number of blocks in the y direction.
         * @override of the abstract function in VolumeTraverser.
         * @return
         */
        [[nodiscard]] U16 GetYCount() const override;

        /**
         * Gets the number of blocks in the z direction.
         * @override of the abstract function in VolumeTraverser.
         * @return
         */
        [[nodiscard]] U16 GetZCount() const override;

        bool IsValidLocation(Vec16 const& block);

      protected:
        bool GoToNextBlock();

        const Domain & mLatticeData;
    };

}

#endif // HEMELB_GEOMETRY_BLOCKTRAVERSER_H
