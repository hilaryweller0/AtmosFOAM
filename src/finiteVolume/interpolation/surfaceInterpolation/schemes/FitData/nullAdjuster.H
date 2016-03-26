#ifndef nullAdjuster_H
#define nullAdjuster_H

#include "adjuster.H"

namespace Foam
{
    class nullAdjuster : public adjuster
    {
    public:
        bool adjustWeights();
    };
}

#endif


