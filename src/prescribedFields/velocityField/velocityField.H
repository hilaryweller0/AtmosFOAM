#ifndef VELOCITY_FIELD_H
#define VELOCITY_FIELD_H

#include "surfaceFields.H"
#include "volFields.H"

namespace Foam
{

class velocityField
{
    public:

    //- Runtime type information
    TypeName("velocityField");

    declareRunTimeSelectionTable
    (
        autoPtr,
        velocityField,
        dict,
        (
            const dictionary& dict
        ),
        (dict)
    );

    static autoPtr<velocityField> New(const dictionary& dict);
    void applyTo(surfaceScalarField& phi) const
    {
        applyTo(phi, phi.time().value());
    }
    void applyTo(surfaceScalarField& phi, scalar time) const;
    virtual ~velocityField() {};

    private:
    virtual void applyToInternalField
    (
        surfaceScalarField& phi,
        scalar time 
    ) const = 0;
    virtual void applyToBoundary
    (
        surfaceScalarField& phi,
        const label patchI,
        scalar time
    ) const = 0;
    virtual void project(surfaceScalarField& phi) const = 0;
};

}

#endif


