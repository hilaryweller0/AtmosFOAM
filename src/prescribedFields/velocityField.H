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
    void applyToInternalField(surfaceScalarField& phi) const;
    void applyToBoundary(surfaceScalarField& phi, const label patchI) const;
    scalar faceFlux(const face& f, const fvMesh& mesh, const Time& t) const;
    virtual vector streamfunctionAt(const point& p, const Time& t) const = 0;
};

}

#endif

