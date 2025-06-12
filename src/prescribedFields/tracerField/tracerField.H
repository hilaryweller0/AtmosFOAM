#ifndef TRACER_FIELD_H
#define TRACER_FIELD_H

#include "advectable.H"
#include "surfaceFields.H"
#include "volFields.H"

namespace Foam
{

class tracerField
{
    private:
    void applyToInternalField(volScalarField& T) const;
    void applyToBoundary(volScalarField& T, const label patchI) const;
    void applyToBoundary(surfaceScalarField& T, const label patchI) const;
    virtual bool hasGradient() const;
    virtual scalar tracerAt(const point& p, const Time& t) const = 0;
    virtual vector gradientAt(const point& p, const Time& t) const;
    const advectable& velocityField;
    const bool spherical_;
    const scalar earthRadius_;

    protected:
    tracerField
    (
        const advectable& velocityField,
        const bool spherical__= false,
        const scalar earthRadius__= 0
    );
    
    public:
    
    //- Runtime type information
    TypeName("tracerField");
    
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

    static autoPtr<tracerField> New
    (
        const dictionary& dict,
        const advectable& velocityField
    );
    void applyTo(volScalarField& T) const;
    void applyTo(surfaceScalarField& Tf) const;
    virtual ~tracerField() {};
    
    bool spherical() const { return spherical_; }
    scalar earthRadius() const { return earthRadius_; }
};

}

#endif


