// The FOAM Project // File: polarPatch.C
/*
-------------------------------------------------------------------------------
 =========         | Class Implementation
 \\      /         |
  \\    /          | Name:   polarPatch
   \\  /           | Family: polarPatch
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION
    Defines a surface mesh with regular latitude-longitude vertices over a
    sphere with un-refinements in the longitude direction so as to maintain a
    more constant spacing of meridians
-------------------------------------------------------------------------------
*/

#include "polarPatch.H"
#include "IOdictionary.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// ************************************************************************* //

// return the list of faces
Foam::faceList Foam::polarPatch::calcPolarFaces
(
    const polarPatchData& grid
) const
{

    // calculate the number of faces
    label nfac = grid.polarCell() ?
        2 - (grid.finestLons().size()-1)
      : 2*(grid.lons()[1].size()-1) - (grid.finestLons().size()-1);

    for(label j=1; j < grid.lons().size()-1; j++)
    {
        nfac += grid.lons()[j].size()-1;
    }

    faceList polarFaces(nfac);

    //Info << "Creating " << nfac << " faces in the polarPatch" << endl;

    // how many vertices does each pole cell have (if it exists)
    label nVertPole = grid.lons()[1].size()+2;

    labelList poleFace(nVertPole);
    labelList tri(3);
    labelList quad(4);
    labelList pent(5);

    label ifac = 0;

    if (grid.polarCell())
    {
        // loop through the vertices surrounding the first pole
        poleFace[0] = 1;
        poleFace[1] = 0;
        for(label i = 2; i < nVertPole; i++)
        {
            poleFace[i] = i;
        }
        //Info << "Face at pole " << ifac << " " << poleFace << endl;
        polarFaces[ifac++] = face(poleFace);
    }
    else
    {
        tri[0] = 0;
        tri[2] = 1;
        for(label i = 0; i < grid.lons()[1].size()-1; i++)
        {
            tri[1] = tri[2];
            tri[2] = i+2;
            //Info << "Face " << ifac << " " << tri << endl;
            polarFaces[ifac++] = face(tri);
        }
    }

    label i0 = 0;
    label i1 = grid.polarCell() ? 2 : 1;

    // loop through the lat-lon faces
    for(label j = 1; j < grid.lats().size()-2; j++)
    {
        i0 = i1;
        i1 = i0 + grid.lons()[j].size();
        if(grid.lons()[j].size() == grid.lons()[j+1].size())
        {
            for(label i = 0; i < grid.lons()[j].size()-1; i++)
            {
                quad[0] = i0 + i;
                quad[1] = i0 + i + 1;
                quad[2] = i1 + i + 1;
                quad[3] = i1 + i;
                //Info << "Face " << ifac << " " << quad << endl;
                polarFaces[ifac++] = face(quad);
            }
        }
        else if(grid.lons()[j].size() < grid.lons()[j+1].size())
        {
            for(label i = 0; i < grid.lons()[j].size()-1; i++)
            {
                pent[0] = i0 + i;
                pent[1] = i0 + i + 1;
                pent[2] = i1 + 2*i + 2;
                pent[3] = i1 + 2*i + 1;
                pent[4] = i1 + 2*i;
                //Info << "Face " << ifac << " " << pent << endl;
                polarFaces[ifac++] = face(pent);
            }
        }
        else
        {
            for(label i = 0; i < grid.lons()[j+1].size()-1; i++)
            {
                pent[0] = i0 + 2*i;
                pent[1] = i0 + 2*i + 1;
                pent[2] = i0 + 2*i + 2;
                pent[3] = i1 + i + 1;
                pent[4] = i1 + i;
                //Info << "Face " << ifac << " " << pent << endl;
                polarFaces[ifac++] = face(pent);
            }
        }
    }

    if (grid.polarCell())
    {
        // loop through the vertices surrounding the second pole
        for(label i = 0; i < nVertPole-2; i++)
        {
            poleFace[i] = i + i1;
        }
        poleFace[nVertPole-2] = i1 + nVertPole-1;
        poleFace[nVertPole-1] = i1 + nVertPole-2;
        //Info << "Face at pole " << ifac << " " << poleFace << endl;
        polarFaces[ifac++] = face(poleFace);
    }
    else
    {
        i0 = i1;
        i1 = i0 + grid.lons()[grid.lats().size()-1].size();
        tri[0] = i1;
        tri[2] = i0;
        for(label i = 0; i < grid.lons()[grid.lats().size()-1].size()-1; i++)
        {
            tri[1] = tri[2];
            tri[2] = i0+i+1;
            //Info << "Face " << ifac << " " << tri << endl;
            polarFaces[ifac++] = face(tri);
        }
    }

    return polarFaces;
}


// ************************************************************************* //

// return the vertex locations
Foam::pointField Foam::polarPatch::calcVerts(const polarPatchData& grid) const
{
    // calculate the number of vertices
    label nvert = 2;
    for(label j=1; j < grid.lons().size()-1; j++)
    {
        nvert += grid.lons()[j].size();
    }
    if (grid.polarCell()) nvert += 2;

    pointField verts(nvert);
    label iv = 0;
    
    verts[iv++] = point(0, grid.lats()[0], 0);
    verts[nvert-1] = point(360, grid.lats()[grid.lats().size()-1], 0);
    if (grid.polarCell())
    {
        verts[iv++] = point(360, grid.lats()[0], 0);
        verts[nvert-2] = point(0, grid.lats()[grid.lats().size()-1], 0);
    }
    
    for(label j = 1; j < grid.lats().size()-1; j++)
    {
        for(label i=0; i < grid.lons()[j].size(); i++)
        {
            verts[iv++] = vector(grid.lons()[j][i], grid.lats()[j], 0);
        }
    }

    return verts;
}

// ************************************************************************* //

void Foam::polarPatch::correctDirection()
{
    for(label ifac = 0; ifac < polarPatch::size(); ifac++)
    {
        scalar dir = polarPatch::operator[](ifac).normal(points()) & vector(0,0,1);
        if(dir < SMALL)
        {
            reverse(operator[](ifac));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polarPatch::polarPatch
(
    const label nlat,
    const label nlon,
    const scalar radius,
    const scalarList& urefLats,
    const scalar maxDxLonRatio,
    const bool polarCell
)
:
    polarPatchData(nlat, nlon, maxDxLonRatio, polarCell),
    PrimitivePatch<face, List, pointField>
    (
        calcPolarFaces(*this),
        calcVerts(*this)
    )
{
    correctDirection();
}

Foam::polarPatch::polarPatch
(
    const IOdictionary& earthProperties
)
:
    polarPatchData(earthProperties),
    PrimitivePatch<face, List, pointField>
    (
        calcPolarFaces(*this),
        calcVerts(*this)
    )
{
    //Info << "points = " << points() << endl;
    //Info << "faces = " << *this << endl;

    correctDirection();
}


// ************************************************************************* //
