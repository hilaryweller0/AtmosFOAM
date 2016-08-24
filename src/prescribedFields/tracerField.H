#ifndef TRACER_FIELD_H
#define TRACER_FIELD_H

#include "fvCFD.H"

namespace Foam
{

class tracerField
{
    public:
    declareRunTimeSelectionTable
    (
        autoPtr,
        tracerField,
        dict,
        (
            const dictionary& dict
        ),
        (dict)
    );

    static autoPtr<tracerField> New(const dictionary& dict);
    void applyTo(volScalarField& T) const;
    virtual ~tracerField() {};
    
    private:
    void applyToInternalField(volScalarField& T) const;
    virtual scalar tracerAt(const point& p, const Time& t) const = 0;
};

}

#endif


