/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::RKfield

Description
    Holds a GeometricField at a number of intermediate steps for use by
    Runge-Kutta time integration with a doubleButcherTableau

SourceFiles
    RKfieldI.H
    RKfield.C
    RKfieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef RKfield_H
#define RKfield_H

#include "GeometricField.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes
class butcherTableau;

/*---------------------------------------------------------------------------*\
                         Class RKfield Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class RKfield
:
    public regIOobject,
    public PtrList<GeometricField<Type, PatchField, GeoMesh> >
{
    // Private data
    
        //- Reference to the Butcher tableau
        const butcherTableau& BT_;
        
        //- Keep track of stage in the Butcher tableau
        label iRK_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RKfield(const RKfield<Type, PatchField, GeoMesh>&);

        //- Disallow default bitwise assignment
        void operator=(const RKfield<Type, PatchField, GeoMesh>&);


public:

    //- Runtime type information
    TypeName("RKfield");

    // Constructors

        //- Construct from IOobject, size, one field and Butcher tableau
        RKfield
        (
            const IOobject&,
            const label nSteps,
            const GeometricField<Type, PatchField, GeoMesh>& f,
            const butcherTableau& BT__
        );

        //- Construct from size, one field and Butcher tableau
        RKfield
        (
            const label nSteps,
            GeometricField<Type, PatchField, GeoMesh>& f,
            const butcherTableau& BT__
        );

    // Member Functions

        //- Do nothing to write out
        bool writeData(Ostream&) const {return true;}

        //- Reference to the Butcher tableau
        const butcherTableau& BT() const {return BT_;}
        
        //- Access/set stage in the Butcher tableau
        label iRK() const { return iRK_;}
        void setRKstage(const label i) { iRK_ = i;}

    // Overloaded operator []
    
        //- Overload operator[] to store iRK_
        GeometricField<Type, PatchField, GeoMesh>& operator[](const label i)
        {
            iRK_ = i;
            return PtrList<GeometricField<Type, PatchField, GeoMesh> >::operator[](i);
        }

        //- Overload operator[] without storing iRK_
        const GeometricField<Type, PatchField, GeoMesh>& operator[]
        (
            const label i
        ) const
        {
            return PtrList<GeometricField<Type, PatchField, GeoMesh> >::operator[](i);
        }
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "RKfield.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
