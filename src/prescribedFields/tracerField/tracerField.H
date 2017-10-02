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
    void applyTo(surfaceScalarField& Tf) const;
    virtual ~tracerField() {};

    protected:
    tracerField(const advectable& velocityField);
    
    private:
    void applyToInternalField(volScalarField& T) const;
    void applyToBoundary(volScalarField& T, const label patchI) const;
    virtual bool hasGradient() const;
    virtual scalar tracerAt(const point& p, const Time& t) const = 0;
    virtual vector gradientAt(const point& p, const Time& t) const;
    const advectable& velocityField;
};

}

#endif


