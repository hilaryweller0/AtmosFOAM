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

\*---------------------------------------------------------------------------*/

#include "FieldToPlot.H"
#include "token.H"
#include "IFstream.H"
#include "stringScalar.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FieldToPlot::FieldToPlot()
:
    name_(),
    plotType_(),
    minMaxDel_(),
    colourScale_()
{}


Foam::FieldToPlot::FieldToPlot
(
    const word& name__,
    const plotTypes& plotType__,
    const FixedList<scalar,3> minMaxDel__,
    const string colourScale__,
    const int dir__
)
:
    name_(name__),
    plotType_(plotType__),
    minMaxDel_(minMaxDel__),
    colourScale_(colourScale__),
    vectorDir_(dir__)
{}


Foam::FieldToPlot::FieldToPlot(Istream& is)
:
    name_(),
    plotType_(),
    minMaxDel_()
{
    operator>>(is, *this);
}

Foam::FieldToPlot::FieldToPlot(const FieldToPlot& ftp)
:
    name_(ftp.name_),
    plotType_(ftp.plotType_),
    minMaxDel_(ftp.minMaxDel_),
    colourScale_(ftp.colourScale_),
    vectorDir_(ftp.vectorDir_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FieldToPlot::~FieldToPlot()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word Foam::FieldToPlot::plotTypeWord() const
{
    return plotType_ == FILLED_CONTOURS ? "filledContours"
         : plotType_ == SOLID_CONTOURS  ? "solidContours"
         : plotType_ == DASHED_CONTOURS ? "dashedContours"
         : plotType_ == VECTORS         ? "vectors"
         : plotType_ == VECTOR_END_POINTS ? "vectorEndPoints"
         : plotType_ == VECTOR_CONTOURS ? "vectorContours"
         : plotType_ == RAW_VALUES      ? "rawValues"
         : plotType_ == RAW_FLUXES      ? "rawFluxes"
         : plotType_ == MESH            ? "mesh"
         : plotType_ == MESHPOINTS      ? "meshPoints"
         : plotType_ == MESHCENTRES     ? "meshCentres"
         : plotType_ == MESH_RANGE       ? "meshRange"
         : plotType_ == ADVECTED_CONTOURS? "advectedContours"
         : plotType_ == NUMBERED        ? "numbered"
         : plotType_ == WRITECONTOURS   ? "writeContours"
         : "none";
}

const Foam::string Foam::FieldToPlot::colour() const
{
    return stringScalar(minMaxDel_[0])/stringScalar(minMaxDel_[1])
           /stringScalar(minMaxDel_[2]);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::operator!=(const Foam::FieldToPlot& ftp1, const Foam::FieldToPlot& ftp2)
{
    return ftp1.name() == ftp2.name()
    && ftp1.plotType() == ftp2.plotType()
    && ftp1.minMaxDel() == ftp2.minMaxDel()
    && ftp1.colourScale() == ftp2.colourScale();
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, FieldToPlot& ftp)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    bool errorState = false;
    if (!t.isPunctuation())
    {
        errorState = true;
    }
    else
    {
        if (t.pToken() != token::BEGIN_SQR)
        {
            errorState = true;
        }
    }
    if (errorState)
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, FieldToPlot&)", is)
            << "wrong token type - expected [ found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Read in the stuff
    word plotType;
    is >> ftp.name_ >> plotType;

    if (plotType == "filledContours")
    {
        ftp.plotType_ = FieldToPlot::FILLED_CONTOURS;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.minMaxDel_[2]
           >> ftp.colourScale_;
    }
    else if (plotType == "contours")
    {
        ftp.plotType_ = FieldToPlot::SOLID_CONTOURS;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.minMaxDel_[2]
            >> ftp.colourScale_;
    }
    else if (plotType == "solidContours")
    {
        ftp.plotType_ = FieldToPlot::SOLID_CONTOURS;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.minMaxDel_[2];
    }
    else if (plotType == "dashedContours")
    {
        ftp.plotType_ = FieldToPlot::DASHED_CONTOURS;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.minMaxDel_[2];
    }
    else if (plotType == "vectors")
    {
        ftp.plotType_ = FieldToPlot::VECTORS;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.colourScale_;
    }
    else if (plotType == "vectorEndPoints")
    {
        ftp.plotType_ = FieldToPlot::VECTOR_END_POINTS;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.colourScale_;
    }
    else if (plotType == "vectorContours")
    {
        ftp.plotType_ = FieldToPlot::VECTOR_CONTOURS;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.minMaxDel_[2]
            >> ftp.colourScale_ >> ftp.vectorDir_;
    }
    else if (plotType == "rawValues")
    {
        ftp.plotType_ = FieldToPlot::RAW_VALUES;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.minMaxDel_[2]
           >> ftp.colourScale_;
    }
    else if (plotType == "rawFluxes")
    {
        ftp.plotType_ = FieldToPlot::RAW_FLUXES;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.minMaxDel_[2]
           >> ftp.colourScale_;
    }
    else if (plotType == "mesh")
    {
        ftp.plotType_ = FieldToPlot::MESH;
        is >> ftp.colourScale_;
    }
    else if (plotType == "meshPoints")
    {
        ftp.plotType_ = FieldToPlot::MESHPOINTS;
        is >> ftp.colourScale_;
    }
    else if (plotType == "meshCentres")
    {
        ftp.plotType_ = FieldToPlot::MESHCENTRES;
        is >> ftp.colourScale_;
    }
    else if (plotType == "meshRange")
    {
        ftp.plotType_ = FieldToPlot::MESH_RANGE;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.colourScale_;
    }
    else if (plotType == "advectedContours")
    {
        ftp.plotType_ = FieldToPlot::ADVECTED_CONTOURS;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.minMaxDel_[2]
            >> ftp.colourScale_;
    }
    else if (plotType == "numbered")
    {
        ftp.plotType_ = FieldToPlot::NUMBERED;
        is >> ftp.minMaxDel_[0];
    }
    else if (plotType == "writeContours")
    {
        ftp.plotType_ = FieldToPlot::WRITECONTOURS;
        is >> ftp.minMaxDel_[0] >> ftp.minMaxDel_[1] >> ftp.minMaxDel_[2]
            >> ftp.colourScale_;
    }
    else
    {
        FatalErrorIn("Istream& operator>>(Istream&, FieldToPlot&)")
        << "Second element of FieldToPlot named " << ftp.name()
        << " should be one of filledContours, solidContours, dashedContours, vectors, vectorContours, mesh, meshPoints, meshCentres, advectedContours, writeContours or rawValues but " << plotType << " given" << exit(FatalError);
    }

    // Check state of IOstream
    is.check("Istream& operator>>(Istream&, FieldToPlot&)");

    is >> t;
    errorState = false;
    if (!t.isPunctuation())
    {
        errorState = true;
    }
    else
    {
        if (t.pToken() != token::END_SQR)
        {
            errorState = true;
        }
    }
    if (errorState)
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, FieldToPlot&)", is)
            << "wrong token type - expected ] found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    return is;
}

Foam::Ostream& Foam::operator<<(Ostream& os, const FieldToPlot& ftp)
{
    os << '[' << ftp.name() << ' ' << ftp.plotTypeWord()
        << ' ' << ftp.minMaxDel()[0] << ' ' << ftp.minMaxDel()[1] << ' '
        << ftp.minMaxDel()[2] << ' ' << ftp.colourScale();
    return os;
}

// ************************************************************************* //
