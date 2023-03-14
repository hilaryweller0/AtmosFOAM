// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / U.name() + ".dat";
{
    OFstream os(tmpFile);
    os << "#x    y    " << U.name() << endl;

    const fvPatchVectorField& Uf = U.boundaryField()[patchi];
    Info << U.name() << " goes from " << min(mag(Uf)) << " to " << max(mag(Uf))
        << endl;

    // Set the centre point and valid distance for Orthographic projection
    point p0(0,0,0);
    scalar maxDist2 = 0;
    scalar r = 0;
    if (ortho)
    {
        point pi = mesh.points()[plotPatch[0][0]];
        r = mag(pi);
        p0 = polarPoint(lon0, lat0, r).cartesian();
        maxDist2 = 2*sqr(r);
    }

    // Find the plot bounding box
    //boundBox bb = setRegion(mesh, projection2D, region, 1.1*r);

    // lon and lat directions if needed
    vectorField rHatf = mesh.boundary()[patchi].Cf();
    vectorField lonHatf = vector(0,0,1) ^ rHatf;
    vectorField latHatf = rHatf ^ lonHatf;
    if (projection2D == "latLon")
    {
        rHatf /= mag(rHatf);
        lonHatf /= mag(lonHatf);
        latHatf /= mag(latHatf);
    }

    for(label faci = 0; faci < Uf.size(); faci+= vectorFreq)
    {
        point p = plotPatch.faceCentres()[faci];

        // First check that vector for this face should be plot
        bool plotFace = !ortho || magSqr(p - p0) < maxDist2;

        if (plotFace)
        {
            if (projection2D == "XY")
            {
                os << p.x() << "   " << p.y() << "   "
                    << p.x() + vectorScale*Uf[faci][0] << "   "
                    << p.y() + vectorScale*Uf[faci][1] << '\n';
            }
            else if (projection2D == "YX")
            {
                os << p.y() << "   " << p.x() << "   "
                    << p.y() + vectorScale*Uf[faci][1] << "   "
                    << p.x() + vectorScale*Uf[faci][0] << '\n';
            }
            else if (projection2D == "XZ")
            {
                os << p.x() << "   " << p.z() << "   "
                    << p.x() + vectorScale*Uf[faci][0] << "   "
                    << p.z() + vectorScale*Uf[faci][2] << '\n';
            }
            else if (projection2D == "ZX")
            {
                os << p.z() << "   " << p.x() << "   "
                    << p.z() + vectorScale*Uf[faci][2] << "   "
                    << p.x() + vectorScale*Uf[faci][0] << '\n';
            }
            else if (projection2D == "YZ")
            {
                os << p.y() << "   " << p.z() << "   "
                    << p.y() + vectorScale*Uf[faci][1] << "   "
                    << p.z() + vectorScale*Uf[faci][2] << '\n';
            }
            else if (projection2D == "ZY")
            {
                os << p.z() << "   " << p.y() << "   "
                    << p.z() + vectorScale*Uf[faci][2] << "   "
                    << p.y() + vectorScale*Uf[faci][1] << '\n';
            }
            else if (projection2D == "latLon")
            {
                scalar r = mag(p);
                scalar sinLat = r > VSMALL ? p.z()/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                scalar cosLat = Foam::sqrt(1 - sqr(sinLat));
                scalar cosLon = mag(r*cosLat) > VSMALL ? p.x()/(r*cosLat):0;
                cosLon = max(min(cosLon, 1), -1);
                scalar lon = radToDeg*sign(p.y())*Foam::acos(cosLon);
                os << lon << "   " << lat << "   "
                   << lon + vectorScale*((Uf[faci]) & (lonHatf[faci])) << "   "
                   << lat + vectorScale*((Uf[faci]) & (latHatf[faci])) << '\n';
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

// plot vectors
systemCall = "gmt psxy '" + tmpFile + "' -J -R -h1 -Sv0.5i+s+e+n"
                  + stringScalar(100./vectorScale)
                  + " -W" + fieldsToPlot[ifield].colourScale()
                  + " -G" + fieldsToPlot[ifield].colourScale()
                  + " -O -K >> " + epsFile;
systemVerbose(systemCall);

// remove temp file
systemCall = "rm '" + tmpFile + "' ";
systemVerbose(systemCall);
