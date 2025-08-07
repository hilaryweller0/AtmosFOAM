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
#include "mathematicalConstants.H"
#include "IOdictionary.H"
#include "dimensionedTypes.H"
#include "Field.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// return the lats for regular spacing
Foam::scalarList Foam::polarPatchData::calcLats(const label nlat) const
{
    scalarList lats(nlat);
    
    const scalar dlat = polarCell_ ? constant::mathematical::pi/scalar(nlat)
                                   : constant::mathematical::pi/scalar(nlat-1);

    lats[0] = polarCell_ ? 0.5*(constant::mathematical::pi - dlat)
                         : 0.5*constant::mathematical::pi;

    for(label ilat = 1; ilat < nlat; ilat++)
    {
        lats[ilat] = lats[0] - ilat*dlat;
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
        latCs[i] = asin(0.5*(sin(lats[i]) + sin(lats[i+1])));
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
    scalarList lons(nlon);
    const scalar dlon = 2*constant::mathematical::pi/scalar(nlon);
    for(label ilon = 0; ilon < nlon; ilon++)
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
    const scalarList& urefLats,
    const scalar ratio,
    const label jEq
) const
{
    List<scalarList> lonsU(lats.size());

    lonsU[jEq] = lons;

    // set the longitudes for the latitudes with j < jEq
    scalar dx = 1;
    scalar dxj = 1;

    scalarList lonsj = lons;
    label lonjSize = lons.size();

    label iLatRef = 0;

    for(label j = jEq-1; j >= 1; j--)
    {
        scalar dlat = mag(lats[j] - lats[j-1]);
        dxj = cos(0.5*(lats[j] + lats[j-1]));
        if
        (
            ((dxj*ratio > SMALL && dxj*ratio < dx)
            ||
            (urefLats.size() > iLatRef
                && mag(lats[j]*180./constant::mathematical::pi) + 0.5*dlat
            >= mag(urefLats[iLatRef])))
        )
        {
            if (lonjSize%2 == 0)
            {
                dx *= 0.5;
                lonjSize /= 2;
                Info << "Unrefinement takes place at latitude = "
                    << lats[j]*scalar(180)/constant::mathematical::pi
                    << " new n lons = " << lonjSize
                    << " exact ratio = " << dx/dxj << endl;;
                for(label i=0; i < lonjSize; i++)
                {
                    lonsj[i] = lonsj[2*i];
                }
                iLatRef++;
            }
            else
            {
                Info << "Cannot un-refine longitude at latitude = "
                    << lats[j]*scalar(180)/constant::mathematical::pi
                    << " because lonjSize = " << lonjSize
                    << " which is not divisible by 2" << endl;
            }
        }

        lonsU[j].setSize(lonjSize);
        for(label i=0; i < lonjSize; i++)
        {
            lonsU[j][i] = lonsj[i];
        }
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
                                      + lons[j2][0] + 2*constant::mathematical::pi);
    }

    return lonCs;
}

// ************************************************************************* //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polarPatchData::polarPatchData
(
    const label nlat,
    const label nlon,
    const scalarList& urefLats,
    const scalar maxDxLonRatio,
    const vector& meshRotation,
    const bool polarCell
)
:
    polarCell_(polarCell),
    maxDxLonRatio_(maxDxLonRatio),
    lats_(calcLats(nlat)),
    latCs_(calcLatCs(lats_)),
    jEq_(calcJEq(lats_)),
    lons_(calcLons(lats_, calcLons(nlon), urefLats, maxDxLonRatio, jEq_)),
    lonCs_(calcLonCs(lons_)),
    meshRotation_(meshRotation)
{
    Info << "vertex latitudes = " << lats_.size() << " ( ";
    for(label i = 0; i < lats_.size(); i++)
    {
        Info << lats_[i]*180/constant::mathematical::pi << ' ';
    }
    Info << " )" << endl;
    Info << "finest vertex longitudes = " << lons_.size() << " ( ";
    for(label i = 0; i < lons_.size(); i++)
    {
        Info << lons_[jEq_][i]*180/constant::mathematical::pi << ' ';
    }
    Info << " )" << endl;
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
        constant::mathematical::pi/180.*scalarField(earthProperties.lookup("latitudes")) :
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
            scalarList(earthProperties.lookup("longitudes"))*constant::mathematical::pi/180.
        ) :
        calcLons
        (
            lats_,
            calcLons(readLabel(earthProperties.lookup("nLon"))),
            earthProperties.lookupOrDefault("urefLats", scalarList(0)),
            earthProperties.lookupOrDefault("maxLonRatio", scalar(0)),
            jEq_
        )
    ),
    lonCs_
    (
        earthProperties.found("longitudes")?
        lons_:
        calcLonCs(lons_)
    ),
    meshRotation_
    (
        earthProperties.lookupOrDefault("meshRotation", vector::zero)
        *constant::mathematical::pi/180.
    )
{
//    if (mag(lats_[0] - latCs_[0]) < VSMALL)
//    {
//        lats_.resize(latCs_.size()+1);
//        lats_[0] = latCs_[0] - 0.5*(latCs_[1] - latCs_[0]);
//        const label n = latCs_.size();
//        for(label j = 1; j < n; j++)
//        {
//            lats_[j] = 0.5*(latCs_[j-1] + latCs_[j]);
//        }
//        lats_[n] = latCs_[n-1] + 0.5*(latCs_[n-2] - latCs_[n-1]);
//    }
//    
//    if (mag(lons_[0][0] - lonCs_[0][0]) < VSMALL)
//    {
//        lons_.resize(lats_.size());
//        for(label j = 0; j < lats_.size(); j++)
//        {
//            lons_[j].setSize(lonCs_[0].size());
//            lons_[j][0] = lonCs_[0][0] - 0.5*(lonCs_[0][1] - lonCs_[0][0]);
//            const label m = lonCs_[0].size();
//            for(label i = 1; i < m; i++)
//            {
//                lons_[j][i] = 0.5*(lonCs_[0][i-1] + lonCs_[0][i]);
//            }
//            lons_[j][m] = lonCs_[0][m-1] + 0.5*(lonCs_[0][m-2] - lonCs_[0][m-1]);
//        }
//    }

    Info << "vertex latitudes = " << lats_.size() << " ( ";
    for(label i = 0; i < lats_.size(); i++)
    {
        Info << lats_[i]*180/constant::mathematical::pi << ' ';
    }
    Info << " )" << endl;
    Info << "finest vertex longitudes = " << lons_.size() << " ( ";
    for(label i = 0; i < lons_.size(); i++)
    {
        Info << lons_[jEq_][i]*180/constant::mathematical::pi << ' ';
    }
    Info << " )" << endl;
}


// ************************************************************************* //
