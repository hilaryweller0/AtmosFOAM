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
#include "polarPoint.H"
#include "IOdictionary.H"
#include "dimensionedTypes.H"
#include "plane.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// ************************************************************************* //

// return the list of faces
Foam::faceList Foam::polarPatch::calcPolarFaces
(
    const polarPatchData& grid
) const
{

    // calculate the number of faces
    label nfac = 0;
    if (grid.polarCell())
    {
        nfac = 2- grid.finestLons().size();
    }
    else
    {
        nfac = - grid.finestLons().size();
    }
    for(label j=0; j < grid.lons().size(); j++)
    {
        nfac += grid.lons()[j].size();
    }

    faceList polarFaces(nfac);

    //Info << "Creating " << nfac << " faces in the polarPatch" << endl;

    // how many vertices does each pole cell have (if it exists)
    label nVertPole = grid.lons()[0].size();

    labelList poleFace(nVertPole);
    labelList tri(3);
    labelList quad(4);
    labelList pent(5);

    label ifac = 0;

    if (grid.polarCell())
    {
        // loop through the vertices surrounding the first pole
        for(label i = 0; i < nVertPole; i++)
        {
            poleFace[i] = i;
        }
        //Info << "Face at pole " << ifac << " " << poleFace << endl;
        polarFaces[ifac++] = face(poleFace);
    }
    else
    {
        tri[0] = 0;
        tri[1] = grid.lons()[0].size();
        tri[2] = 1;
        //Info << "Face " << ifac << " " << tri << endl;
        polarFaces[ifac++] = face(tri);
        for(label i = 0; i < grid.lons()[0].size()-1; i++)
        {
            tri[1] = tri[2];
            tri[2] = i+2;
            //Info << "Face " << ifac << " " << tri << endl;
            polarFaces[ifac++] = face(tri);
        }
    }

    label i0 = 0;
    label i1 = grid.polarCell() ? 0 : 1;

    // loop through the lat-lon faces
    const label j0 = grid.polarCell() ? 0 : 1;
    for(label j = j0; j < grid.lats().size()-j0-1; j++)
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
            quad[0] = i1-1;
            quad[1] = i0;
            quad[2] = i1;
            quad[3] = i1 + grid.lons()[j].size() - 1;
            //Info << "Face " << ifac << " " << quad << endl;
            polarFaces[ifac++] = face(quad);
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
            if(grid.lons()[j+1].size()/2 == grid.lons()[j].size())
            {
                pent[0] = i1-1;
                pent[1] = i0;
                pent[2] = i1;
                pent[3] = i1 + grid.lons()[j+1].size() - 1;
                pent[4] = i1 + grid.lons()[j+1].size() - 2;
                //Info << "Face " << ifac << " " << pent << endl;
                polarFaces[ifac++] = face(pent);
            }
            else
            {
                FatalErrorIn("Foam::polarPatch::calcPolarFaces")
                    << " can only have quads and pents in polarPatch"
                    << exit(FatalError);
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
            if(grid.lons()[j].size()/2 == grid.lons()[j+1].size())
            {
                pent[0] = i1-2;
                pent[1] = i1-1;
                pent[2] = i0;
                pent[3] = i1;
                pent[4] = i1 + grid.lons()[j+1].size() - 1;
                //Info << "Face " << ifac << " " << pent << endl;
                polarFaces[ifac++] = face(pent);
            }
            else
            {
                FatalErrorIn("Foam::polarPatch::calcPolarFaces")
                    << " can only have quads and pents in polarPatch"
                    << exit(FatalError);
            }
        }
    }

    if (grid.polarCell())
    {
        // loop through the vertices surrounding the second pole
        for(label i = 0; i < nVertPole; i++)
        {
            poleFace[i] = i + i1;
        }
        //Info << "Face at pole " << ifac << " " << poleFace << endl;
        polarFaces[ifac++] = face(poleFace);
    }
    else
    {
        i0 = i1;
        i1 = i0 + grid.lons()[grid.lats().size()-1].size();
        tri[0] = i1;
        tri[1] = i1-1;
        tri[2] = i0;
        //Info << "Face " << ifac << " " << tri << endl;
        polarFaces[ifac++] = face(tri);
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
Foam::pointField Foam::polarPatch::calcVerts
(
    const polarPatchData& grid,
    const scalar r
) const
{
    // calculate the number of vertices
    label nvert = 0;
    if (!grid.polarCell()) nvert = 2*(1 - grid.lons()[0].size());
    for(label j=0; j < grid.lons().size(); j++)
    {
        nvert += grid.lons()[j].size();
    }

    pointField verts(nvert);
    const vector& a = grid.meshRotation();

    // rotation tensor:
    tensor rot
    (
        1, 0,         0,
        0, cos(a.x()),-sin(a.x()),
        0, sin(a.x()), cos(a.x())
    );

    rot = rot & tensor
    (
        cos(a.y()), 0,-sin(a.y()),
        0,          1, 0,
        sin(a.y()), 0, cos(a.y())
    );
    rot = rot & tensor
    (
        cos(a.z()), -sin(a.z()), 0,
        sin(a.z()),  cos(a.z()), 0,
        0,  0 , 1
    );
    
    Info << "Rotation vector " << a << endl;
    Info << "Rotation tensor " << rot << endl;

    label iv = 0;
    
    if (!grid.polarCell())
    {
        verts[iv] = polarPoint(0, grid.lats()[0], r).cartesian();
        verts[iv] = rot & (verts[iv]);
        iv++;
    }
    
    const label j0 = grid.polarCell() ? 0 : 1;
    for(label j = j0; j < grid.lats().size()-j0; j++)
    {
        for(label i=0; i < grid.lons()[j].size(); i++)
        {
            verts[iv] = polarPoint
                        (
                            grid.lons()[j][i],
                            grid.lats()[j],
                            r
                        ).cartesian();
            verts[iv] = rot & (verts[iv]);
            iv++;
        }
    }

    if (!grid.polarCell())
    {
        verts[iv]=polarPoint(0,grid.lats()[grid.lats().size()-1],r).cartesian();
        verts[iv] = rot & (verts[iv]);
        iv++;
    }

    return verts;
}

// ************************************************************************* //

void Foam::polarPatch::correctDirection()
{
    for(label ifac = 0; ifac < polarPatch::size(); ifac++)
    {
        scalar dir = polarPatch::operator[](ifac).normal(points())
                   & polarPatch::operator[](ifac).centre(points());
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
    const vector& meshRotation,
    const bool polarCell
)
:
    polarPatchData(nlat, nlon, urefLats, maxDxLonRatio, meshRotation,polarCell),
    PrimitivePatch<face, List, pointField>
    (
        calcPolarFaces(*this),
        calcVerts(*this, radius)
    ),
    radius_(radius)
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
        calcVerts(*this, dimensionedScalar(earthProperties.lookup("earthRadius")).value())
    ),
    radius_(dimensionedScalar(earthProperties.lookup("earthRadius")).value())
{
    correctDirection();
}


Foam::vectorField Foam::polarPatch::voronoiPoints() const
{
    vectorField vPts(size(), vector::zero);
    
    forAll(vPts, cellI)
    {
        const face& f = (*this)[cellI];
        scalar n = 0;
        forAll(f, i)
        {
            vPts[cellI] += points()[f[i]];
            n++;
            if (!polarCell() && (f[i] == 0 || f[i] == nPoints()-1))
            {
                vPts[cellI] += points()[f[i]];
                n++;
            }
        }
        vPts[cellI] /= n;
    }
    return vPts;
}

Foam::vectorField Foam::polarPatch::edgePoints(const vectorField& vPts) const
{
    vectorField ePts(nEdges(), vector::zero);
    
    forAll(ePts, ei)
    {
        const edge& e = edges()[ei];
        const point& p0 = localPoints()[e[0]];
        const point& p1 = localPoints()[e[1]];
        const point& v0 = vPts[edgeFaces()[ei][0]];
        const point& v1 = vPts[edgeFaces()[ei][1]];
        
        const plane edgePlane(vector::zero, p0, p1);
        const plane vPlane(vector::zero, v0, v1);
        const plane::ray r = edgePlane.planeIntersect(vPlane);
        
        ePts[ei] = r.dir()/mag(r.dir())*mag(p0)*sign(r.dir() & p0);
    }
    return ePts;
}

// ************************************************************************* //

