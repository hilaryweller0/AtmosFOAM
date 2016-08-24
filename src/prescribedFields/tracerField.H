#ifndef TRACER_FIELD_H
#define TRACER_FIELD_H

#include "fvCFD.H"
#include "advectable.H"

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
            const dictionary& dict,
            const advectable& velocityField
        ),
        (dict, velocityField)
    );

    static autoPtr<tracerField> New(const dictionary& dict, const advectable& velocityField);
    void applyTo(volScalarField& T) const;
    virtual ~tracerField() {};

    protected:
    tracerField(const advectable& velocityField);
    
    private:
    void applyToInternalField(volScalarField& T) const;
    virtual scalar tracerAt(const point& p, const Time& t) const = 0;
    const advectable& velocityField;
};

}

#endif


