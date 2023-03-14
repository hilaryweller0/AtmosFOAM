/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::epsilonCmuWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulence dissipation wall constraint
    for low- and high-Reynolds number turbulence models.
    This version must be used with stratifiedBoundaryFoam and with
    nutkAtmRoughCmuWallFunction as it looks up Cmu from the database.

    The condition can be applied to wall boundaries for which it
    - calculates \c epsilon and \c G
    - specifies the near-wall epsilon value

    where

    \vartable
        epsilon | turbulence dissipation field
        G       | turbulence generation field
    \endvartable

    The model switches between laminar and turbulent functions based on the
    laminar-to-turbulent y+ value derived from the kappa and E specified in the
    corresponding nutWallFunction.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            epsilonCmuWallFunction;
    }
    \endverbatim

See also
    Foam::fixedInternalValueFvPatchField
    Foam::omegaWallFunctionFvPatchScalarField

SourceFiles
    epsilonCmuWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonCmuWallFunctionFvPatchScalarField_H
#define epsilonCmuWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class momentumTransportModel;

/*---------------------------------------------------------------------------*\
           Class epsilonCmuWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class epsilonCmuWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

        //- Tolerance used in weighted calculations
        static scalar tolerance_;

        //- Local copy of turbulence G field
        scalarField G_;

        //- Local copy of turbulence epsilon field
        scalarField epsilon_;

        //- Initialised flag
        bool initialised_;

        //- Master patch ID
        label master_;

        //- List of averaging corner weights
        List<List<scalar>> cornerWeights_;


    // Protected Member Functions

        //- Set the master patch - master is responsible for updating all
        //  wall function patches
        virtual void setMaster();

        //- Create the averaging weights for cells which are bounded by
        //  multiple wall function faces
        virtual void createAveragingWeights();

        //- Helper function to return non-const access to an epsilon patch
        virtual epsilonCmuWallFunctionFvPatchScalarField& epsilonPatch
        (
            const label patchi
        );

        //- Main driver to calculate the turbulence fields
        virtual void calculateTurbulenceFields
        (
            const momentumTransportModel& turbModel,
            scalarField& G0,
            scalarField& epsilon0
        );

        //- Calculate the epsilon and G
        virtual void calculate
        (
            const momentumTransportModel& turbModel,
            const List<scalar>& cornerWeights,
            const fvPatch& patch,
            scalarField& G,
            scalarField& epsilon
        );

        //- Return non-const access to the master patch ID
        virtual label& master()
        {
            return master_;
        }


public:

    //- Runtime type information
    TypeName("epsilonCmuWallFunction");


    // Constructors

        //- Construct from patch and internal field
        epsilonCmuWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        epsilonCmuWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  epsilonCmuWallFunctionFvPatchScalarField
        //  onto a new patch
        epsilonCmuWallFunctionFvPatchScalarField
        (
            const epsilonCmuWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        epsilonCmuWallFunctionFvPatchScalarField
        (
            const epsilonCmuWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        epsilonCmuWallFunctionFvPatchScalarField
        (
            const epsilonCmuWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new epsilonCmuWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~epsilonCmuWallFunctionFvPatchScalarField()
    {}


    // Member Functions

        // Access

            //- Return non-const access to the master's G field
            scalarField& G(bool init = false);

            //- Return non-const access to the master's epsilon field
            scalarField& epsilon(bool init = false);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
