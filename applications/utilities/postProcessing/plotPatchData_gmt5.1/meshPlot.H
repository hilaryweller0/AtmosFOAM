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
    point pi = mesh.points()[0];

    scalar r = mag(pi);
    if (ortho)
    {
        p0 = polarPoint(lon0, lat0, r).cartesian();
        maxDist2 = 2*sqr(r);
       //Info << "p0 = " << p0 << " maxDist = " << Foam::sqrt(maxDist2) <<endl;
    }
    // Find the plot bounding box
    //boundBox bb = setRegion(mesh, projection2D, region, r);

    for(label faci = 0; faci < plotPatch.size(); faci++)
    {
        face f = plotPatch[faci];

        // First check that all vertices of this face should be plot
        bool plotFace = true;
        if (ortho)
        {
            for(label ip = 0; ip < f.size() && plotFace; ip++)
            {
                point p = mesh.points()[f[ip]];
                if (magSqr(p - p0) >= maxDist2) plotFace = false;
            }
        }

        if (plotFace)
        {
            os << ">\n";
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
                else if (projection2D == "latLon")
                {
                    scalar r = mag(p);
                    scalar sinLat = r > VSMALL ? p.z()/r : 0;
                    sinLat = max(min(sinLat, 1), -1);
                    scalar lat = radToDeg*Foam::asin(sinLat);
                    scalar cosLat = Foam::sqrt(1 - sqr(sinLat));
                    scalar cosLon = mag(r*cosLat) > VSMALL ?
                    p.x()/(r*cosLat):0;
                    cosLon = max(min(cosLon, 1), -1);
                    scalar lon = radToDeg*sign(p.y())*Foam::acos(cosLon);
                    os << setprecision(12) << lon << "   " << lat << '\n';
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
    }
    os << endl;
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
