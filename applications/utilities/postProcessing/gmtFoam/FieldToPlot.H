/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    FieldToPlot

Description

SourceFiles
    FieldToPlot.C

\*---------------------------------------------------------------------------*/

#ifndef FieldToPlot_H
#define FieldToPlot_H

#include "word.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class FieldToPlot;
bool operator!=(const FieldToPlot&, const FieldToPlot&);
Ostream& operator<<(Ostream&, const FieldToPlot&);
Istream& operator>>(Istream&, FieldToPlot&);

/*---------------------------------------------------------------------------*\
                           Class FieldToPlot Declaration
\*---------------------------------------------------------------------------*/

class FieldToPlot
{
public:

    enum plotTypes
    {
        FILLED_CONTOURS,
        SOLID_CONTOURS,
        DASHED_CONTOURS,
        VECTORS,
        VECTOR_END_POINTS,
        VECTOR_CONTOURS,
        RAW_VALUES,
        RAW_FLUXES,
        MESH,
        MESHPOINTS,
        MESHCENTRES,
        MESH_RANGE,
        ADVECTED_CONTOURS,
        NUMBERED,
        WRITECONTOURS,
        DEFAULT
    };

    // Private data

    word name_;
    plotTypes plotType_;
    FixedList<scalar,3> minMaxDel_;
    string colourScale_;
    int vectorDir_;

    // Private Member Functions

public:

    // Constructors

        //- Construct null
        FieldToPlot();

        //- Construct from components
        FieldToPlot
        (
            const word& name__,
            const plotTypes& plotType__,
            const FixedList<scalar,3> minMaxDel__,
            const string colourScale__ = "",
            const int dir__ = 0
        );

        //- Construct from Istream
        FieldToPlot(Istream&);

        //- Construct as copy
        FieldToPlot(const FieldToPlot&);


    // Destructor

        ~FieldToPlot();


    // Member Functions

        // Access

        const word& name() const {return name_;}
        plotTypes plotType() const {return plotType_;}
        const word plotTypeWord() const;
        const FixedList<scalar,3>& minMaxDel() const {return minMaxDel_;}
        const string& colourScale() const {return colourScale_;}
        scalar min() const { return minMaxDel_[0]; }
        scalar max() const { return minMaxDel_[1]; }
        scalar delta() const { return minMaxDel_[2]; }
        scalar vectorScale() const { return minMaxDel_[0]; }
        label  vectorFreq()  const { return label(minMaxDel_[1]); }
        label fontSize() const { return label(minMaxDel_[0]); }
        const string colour() const;
        int vectorDir() const { return vectorDir_; }

    // Friend Operators

        friend bool operator!=(const FieldToPlot& ftp1,const FieldToPlot& ftp2);

    // IOstream Operators

        friend Istream& operator>>(Istream&, FieldToPlot&);
        friend Ostream& operator<<(Ostream&, const FieldToPlot&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
