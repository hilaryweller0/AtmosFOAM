// The FOAM Project // File: polarPatchData.C
/*
-------------------------------------------------------------------------------
 =========         | Class Implementation
 \\      /         |
  \\    /          | Name:   polarPatchData
   \\  /           | Family: polarPatch
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION
    Defines the regular latitude-longitude locations over a
    sphere with un-refinements in the longitude direction so as to maintain a
    more constant spacing of meridians

AUTHOR
    Hilary Spencer.

-------------------------------------------------------------------------------
*/

#include "polarPatchData.H"
#include "IOdictionary.H"
#include "dimensionedTypes.H"
#include "Field.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// return the lats for regular spacing
Foam::scalarList Foam::polarPatchData::calcLats(const label nlat) const
{
    scalarList lats = polarCell_ ? scalarList(nlat+2) : scalarList(nlat+1);
    
    const scalar dlat = 180./scalar(nlat);

    lats[0] = 90;
    lats[1] = polarCell_ ? 90 - 0.5*dlat : 90 - dlat;
    lats[lats.size()-1] = -90;
    lats[lats.size()-2] = polarCell_ ? -90 + 0.5*dlat : -90 + dlat;

    for(label ilat = 2; ilat < lats.size()-2; ilat++)
    {
        lats[ilat] = lats[1] - (ilat-1)*dlat;
    }

    return lats;
}

// ************************************************************************* //

// return the lat cell centres
Foam::scalarList Foam::polarPatchData::calcLatCs(const scalarList& lats) const
{
    scalarList latCs(lats.size()-1);
    for(label i=0; i < latCs.size(); i++)
    {
        latCs[i] = 0.5*(lats[i] + lats[i+1]);
    }
    return latCs;
}

// ************************************************************************* //

// calculate the label for the latitude closest to the equator
Foam::label Foam::polarPatchData::calcJEq(const scalarList& lats) const
{
    label jEq = 0;
    scalar latEq = mag(lats[jEq]);

    for(label j=1; j < lats.size(); j++)
    {
        if(mag(lats[j]) < latEq)
        {
            latEq = mag(lats_[j]);
            jEq = j;
        }
    }
    return jEq;
}

// ************************************************************************* //

// return the lons for regular spacing
Foam::scalarList Foam::polarPatchData::calcLons(const label nlon) const
{
    scalarList lons(nlon+1);
    const scalar dlon = 360./scalar(nlon);
    for(label ilon = 0; ilon <= nlon; ilon++)
    {
        lons[ilon] = ilon*dlon;
    }
    return lons;
}

// ************************************************************************* //

// set the longitudes with unrefinement
Foam::List<Foam::scalarList> Foam::polarPatchData::calcLons
(
    const scalarList& lats,
    const scalarList& lons,
    const scalar maxDxLonRatio,
    const label jEq
) const
{
    List<scalarList> lonsU(lats.size());

    lonsU[jEq] = lons;

    // set the longitudes for the latitudes with j < jEq
    scalarList lonsj = lons;
    label lonjSize = lons.size();
    scalar ratio = 1;
    label jLatRef = jEq;

    for(label j = jEq-1; j >= 1; j--)
    {
        scalar cosLat = cos(lats[j]*constant::mathematical::pi/180.);
        if (ratio/cosLat > maxDxLonRatio && j != 1 && j+1 != jLatRef)
        {
            if ((lonjSize-1)%2 == 0)
            {
                ratio *= 0.5;
                jLatRef = j;
                lonjSize = (lonjSize-1)/2 + 1;
                Info << "Unrefinement takes place at latitude = "
                    << lats[j] << " new n lons = " << lonjSize << endl;
                for(label i=0; i < lonjSize; i++)
                {
                    lonsj[i] = lonsj[2*i];
                }
                Info << "New longitude locations = " << lonjSize << "\n(";
                for(label i = 0; i < lonjSize; i++) Info << lonsj[i] << " ";
                Info << ")" << endl;
            }
            else
            {
                Info << "Cannot un-refine longitude at latitude = "
                    << lats[j] << " because lonjSize = " << lonjSize
                    << " which is divisible by 2" << endl;
            }
        }

        lonsU[j].setSize(lonjSize);
        for(label i = 0; i < lonjSize; i++) lonsU[j][i] = lonsj[i];
    }
    lonsU[0] = lonsU[1];

    // set the longitudes for the latitudes with j > jEq
    for(label j = jEq+1; j < lats.size(); j++)
    {
        lonsU[j] = lonsU[lats.size()-1-j];
    }

    return lonsU;
}

// ************************************************************************* //

// return the lon cell centres
Foam::List<Foam::scalarList> Foam::polarPatchData::calcLonCs
(
    const List<scalarList>& lons
) const
{
    List<scalarList> lonCs(lons.size() - 1);

    for(label j = 0; j < lonCs.size(); j++)
    {
        label j2 = j;
        if(lons[j+1].size() < lons[j].size())
        {
            j2 = j+1;
        }
        lonCs[j].setSize(lons[j2].size());

        for(label i = 0; i < lonCs[j].size()-1; i++)
        {
            lonCs[j][i] = 0.5*(lons[j2][i] + lons[j2][i+1]);
        }
        lonCs[j][lonCs[j].size()-1] = 0.5*(lons[j2][lonCs[j].size()-1]
                                      + lons[j2][0] + 360.);
    }

    return lonCs;
}

// ************************************************************************* //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polarPatchData::polarPatchData
(
    const label nlat,
    const label nlon,
    const scalar maxDxLonRatio,
    const bool polarCell
)
:
    polarCell_(polarCell),
    maxDxLonRatio_(maxDxLonRatio),
    lats_(calcLats(nlat)),
    latCs_(calcLatCs(lats_)),
    jEq_(calcJEq(lats_)),
    lons_(calcLons(lats_, calcLons(nlon), maxDxLonRatio, jEq_)),
    lonCs_(calcLonCs(lons_))
{
    Info << "vertex latitudes = " << lats_.size() << " ( ";
    for(label i = 0; i < lats_.size(); i++)
    {
        Info << lats_[i] << ' ';
    }
    Info << " )" << endl;
    Info << "finest vertex longitudes = " << finestLons().size() << " ( ";
    for(label i = 0; i < finestLons().size(); i++)
    {
        Info << finestLons()[i] << ' ';
    }
    Info << " )" << endl;
    
//    Info << "All vertex longitudes" << endl;
//    for(label j = 0; j < lons_.size(); j++)
//    {
//        Info << lons_[j].size() << " (";
//        for(label i = 0; i < lons_[j].size(); j++)
//        {
//            Info << lons_[j][j] << " ";
//        }
//        Info << ")" << endl;
//    }
}


Foam::polarPatchData::polarPatchData
(
    const IOdictionary& earthProperties
)
:
    polarCell_(earthProperties.lookupOrDefault("polarCell", true)),
    maxDxLonRatio_(earthProperties.lookupOrDefault("maxLonRatio", scalar(0))),
    lats_
    (
        earthProperties.found("latitudes")?
        scalarField(earthProperties.lookup("latitudes")) :
        calcLats(readLabel(earthProperties.lookup("nLat")))
    ),
    latCs_
    (
        earthProperties.found("latitudes")?
        lats_ :
        calcLatCs(lats_)
    ),
    jEq_(calcJEq(lats_)),
    lons_
    (
        earthProperties.found("longitudes")?
        List<scalarList>
        (
            lats_.size(),
            scalarList(earthProperties.lookup("longitudes"))
        ) :
        calcLons
        (
            lats_,
            calcLons(readLabel(earthProperties.lookup("nLon"))),
            earthProperties.lookupOrDefault("maxLonRatio", scalar(0)),
            jEq_
        )
    ),
    lonCs_
    (
        earthProperties.found("longitudes")?
        lons_:
        calcLonCs(lons_)
    )
{
    Info << "vertex latitudes = " << lats_.size() << " ( ";
    for(label i = 0; i < lats_.size(); i++)
    {
        Info << lats_[i] << ' ';
    }
    Info << " )" << endl;
    Info << "finest vertex longitudes = " << finestLons().size() << " ( ";
    for(label i = 0; i < finestLons().size(); i++)
    {
        Info << finestLons()[i] << ' ';
    }
    Info << " )" << endl;

    //Info << "All vertex longitudes = " << lons_ << endl;
}


// ************************************************************************* //
