#ifndef VELOCITY_FIELD_H
#define VELOCITY_FIELD_H

#include "fvCFD.H"

namespace Foam
{

class velocityField
{
    public:
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
    void applyTo(surfaceScalarField& phi) const;
    virtual ~velocityField() {};

    private:
    virtual void applyToInternalField(surfaceScalarField& phi) const = 0;
    virtual void applyToBoundary(surfaceScalarField& phi, const label patchI) const = 0;
    virtual void project(surfaceScalarField& phi) const = 0;
};

}

#endif


