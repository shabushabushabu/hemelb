// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/writers/XdrWriter.h"

namespace hemelb::io
{

        // No field/record separators in XDR files
        void XdrWriter::writeFieldSeparator()
        {
        }
        void XdrWriter::writeRecordSeparator()
        {
        }

}
