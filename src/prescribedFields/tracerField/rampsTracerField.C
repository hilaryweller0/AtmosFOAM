#include "rampsTracerField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
using namespace constant::mathematical;

defineTypeNameAndDebug(rampsTracerField, 0);
addToRunTimeSelectionTable(tracerField, rampsTracerField, dict);

rampsTracerField::rampsTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    nRamps_(readLabel(dict.lookup("nRamps"))),
    centres_(nRamps_),
    lengths_(nRamps_),
    maxTracers_(nRamps_)
{
    for(label i = 0; i < nRamps_; i++)
    {
        centres_[i] = dict.lookup(word("rampCentre"+std::to_string(i)));
        lengths_[i] = dict.lookup(word("rampLength"+std::to_string(i)));
        maxTracers_[i]
            = readScalar(dict.lookup(word("rampMax"+std::to_string(i))));
    }
}

scalar rampsTracerField::tracerAt(const point& p, const Time& t) const
{
    scalar tracer = 0;
    
    // Add contribution from each ramp
    for(label ir = 0; ir < nRamps_; ir++)
    {
        // Vector to ramp centre
        vector d = p - centres_[ir];
        scalar dist = mag(d & lengths_[ir])/(lengths_[ir] & lengths_[ir]);
        if (dist < 1)
        {
            tracer += maxTracers_[ir]*sqr(Foam::cos(0.5*pi*dist));
        }
    }

    return tracer;
}

