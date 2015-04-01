{
// temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / "mesh.dat";
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

    for(label celli = 0; celli < mesh.nCells(); celli++)
    {
        // First check that all vertices of this cell should be plot
        bool plotCell = true;

        for(label faci = 0; faci < mesh.cells()[celli].size()&&plotCell;faci++)
        {
            face f = mesh.faces()[mesh.cells()[celli][faci]];

            // Check that all vertices of this face should be plot
            if (ortho)
            {
                for(label ip = 0; ip < f.size() && plotCell; ip++)
                {
                    point p = mesh.points()[f[ip]];
                    if (magSqr(p - p0) >= maxDist2) plotCell = false;
                }
            }
            
            // Check that the cell centre is within the ranges
            const point& p = mesh.C()[celli];
            if (plotCell && rMax > rMin)
            {
                const scalar r = mag(p);
                if (r < rMin || r > rMax) plotCell = false;
            }
            if (plotCell && lonMax > lonMin)
            {
                const scalar lon = radToDeg*Foam::atan2((p&ydir),(p&xdir));
                if (lon < lonMin || lon > lonMax) plotCell = false;
            }
            if (plotCell && latMax > latMin)
            {
                const scalar lat = radToDeg*Foam::asin(p&axis/mag(p));
                if (lat < latMin || lat > latMax) plotCell = false;
            }
        }

        if (plotCell)
        {
            for(label faci = 0; faci < mesh.cells()[celli].size(); faci++)
            {
                os << ">\n";
                const label faceI = mesh.cells()[celli][faci];
                face f = mesh.faces()[faceI];

//                if (faceI < mesh.nInternalFaces())
//                {
//                    const point& p = mesh.Cf()[faceI];
//                    OStringStream ostrS;
//                    ostrS << p.x() << ' ' << p.y() << " 12 0 4 CM " << faceI;
//                    systemCall = "echo " + ostrS.str()
//                                      + " | pstext  -J -R -h0 -K -O >>"
//                                      + epsFile;
//                    systemVerbose(systemCall);
//                }

                forAll (f, ip)
                {
                    point p = mesh.points()[f[ip]];

                    if (projection2D == "XY")
                    {
                        os << p.x() << "   " << p.y() << '\n';
                    }
                    else if (projection2D == "YX")
                    {
                        os << p.y() << "   " << p.x() << '\n';
                    }
                    else if (projection2D == "XZ")
                    {
                        os << p.x() << "   " << p.z() << '\n';
                    }
                    else if (projection2D == "ZX")
                    {
                        os << p.z() << "   " << p.x() << '\n';
                    }
                    else if (projection2D == "YZ")
                    {
                        os << p.y() << "   " << p.z() << '\n';
                    }
                    else if (projection2D == "ZY")
                    {
                        os << p.z() << "   " << p.y() << '\n';
                    }
                    else if (projection2D == "latz")
                    {
                        scalar r = mag(p);
                        scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                        sinLat = max(min(sinLat, 1), -1);
                        scalar lat = radToDeg*Foam::asin(sinLat);
                        os << lat << "   " << r << '\n';
                    }
                    else if (projection2D == "lonz")
                    {
                        scalar r = mag(p);
                        scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
                        os << lon << "   " << r << '\n';
                    }
                    else if (projection2D == "latLon")
                    {
                        scalar r = mag(p);
                        scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                        sinLat = max(min(sinLat, 1), -1);
                        scalar lat = radToDeg*Foam::asin(sinLat);
                        scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
                        os << lon << "   " << lat << '\n';
                    }
                    else
                    {
                        FatalErrorIn("gmtFoam")
                    << "projection2D should be one of XY, XZ, YZ or latLon but "
                        << projection2D << " given"
                        << exit(FatalError);
                    }
                }
            }
        }
    }
}

// plot mesh
systemCall = "psxy '" + tmpFile + "' -W"
                  + fieldsToPlot[ifield].colourScale()
                           + " -J -R -h1 -L -K -O >>" + epsFile;
systemVerbose(systemCall);

// delete temporary files
systemCall = "rm '" + tmpFile + "' ";
systemVerbose(systemCall);
}
