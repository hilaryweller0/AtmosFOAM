{
// temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / "boundary.dat";
{
    OFstream os(tmpFile);
    os << "#x    y" << endl;

    // Set the centre point and valid distance for Orthographic projection
    point p0(0,0,0);
    scalar maxDist2 = 0;
    point pi = mesh.C()[0];

    scalar r = mag(pi);
    if (ortho)
    {
//        p0 = polarPoint(lon0, lat0, r).cartesian(xdir, ydir, axis);
        p0 = polarPoint(lon0, lat0, r).cartesian();
        maxDist2 = 2.5*sqr(r);
    }

    // Loop over the boundaries and then over the edges of each boundary
    for(label patchI = 0; patchI < mesh.boundaryMesh().size(); patchI++)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];
        if (!isA<emptyPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
        forAll(patch.edgeLoops(), loopi)
        {
            os << "\n>\n";
            forAll(patch.edgeLoops()[loopi], ip)
            {
                const point& pt = patch.localPoints()[patch.edgeLoops()[loopi], ip];
            
                // First check that this vertex should be plot
                bool plotCell = true;

                // Check that all vertices of this face should be plot
                if (ortho)
                {
                    if (magSqr(pt - p0) >= maxDist2) plotCell = false;
                }
                
                // Check that this point is within the ranges
                if (plotCell && rMax > rMin)
                {
                    const scalar r = mag(pt);
                    if (r < rMin || r > rMax) plotCell = false;
                }
                if (plotCell && lonMax > lonMin)
                {
                    const scalar lon = radToDeg*Foam::atan2((pt&ydir),(pt&xdir));
                    if (lon < lonMin || lon > lonMax) plotCell = false;
                }
                if (plotCell && latMax > latMin)
                {
                    const scalar lat = radToDeg*Foam::asin(pt&axis/mag(pt));
                    if (lat < latMin || lat > latMax) plotCell = false;
                }

                if (plotCell)
                {
                    if (projection2D == "XY")
                    {
                        os << pt.x() << "   " << pt.y() << '\n';
                    }
                    else if (projection2D == "YX")
                    {
                        os << pt.y() << "   " << pt.x() << '\n';
                    }
                    else if (projection2D == "XZ")
                    {
                        os << pt.x() << "   " << pt.z() << '\n';
                    }
                    else if (projection2D == "ZX")
                    {
                        os << pt.z() << "   " << pt.x() << '\n';
                    }
                    else if (projection2D == "YZ")
                    {
                        os << pt.y() << "   " << pt.z() << '\n';
                    }
                    else if (projection2D == "ZY")
                    {
                        os << pt.z() << "   " << pt.y() << '\n';
                    }
                    else if (projection2D == "latz")
                    {
                        scalar r = mag(pt);
                        scalar sinLat = r > VSMALL ? (pt&axis)/r : 0;
                        sinLat = max(min(sinLat, 1), -1);
                        scalar lat = radToDeg*Foam::asin(sinLat);
                        os << lat << "   " << r << '\n';
                    }
                    else if (projection2D == "lonz")
                    {
                        scalar r = mag(pt);
                        scalar lon = radToDeg*Foam::atan2((pt&ydir), (pt&xdir));
                        os << lon << "   " << r << '\n';
                    }
                    else if (projection2D == "latLon")
                    {
                        scalar r = mag(pt);
                        scalar sinLat = r > VSMALL ? (pt&axis)/r : 0;
                        sinLat = max(min(sinLat, 1), -1);
                        scalar lat = radToDeg*Foam::asin(sinLat);
                        scalar lon = radToDeg*Foam::atan2((pt&ydir), (pt&xdir));
                        os << lon << "   " << lat << '\n';
                    }
                    else
                    {
                        FatalErrorIn("plotPatchData")
                    << "projection2D should be one of XY, XZ, YZ or latLon but "
                        << projection2D << " given"
                        << exit(FatalError);
                    }
                }
            }
        }}
    }
}

// plot boundary
systemCall = "psxy '" + tmpFile + "' -Sp2p"
                      + " -J -R -h1 -K -O >>" + epsFile;
systemVerbose(systemCall);

// delete temporary files
systemCall = "rm '" + tmpFile + "' ";
systemVerbose(systemCall);
}
