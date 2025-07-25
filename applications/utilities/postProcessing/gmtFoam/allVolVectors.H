// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.name() / U.name() + ".dat";
scalar maxMagU = 0;
{
    OFstream os(tmpFile);
    os << "#x    y    " << U.name() << endl;

    const DimensionedField<vector,volMesh>& Uf = U.internalField();
    maxMagU = max(mag(Uf)).value();
    Info << U.name() << " goes from " << min(mag(Uf)).value() << " to "
         << maxMagU << endl;

    if (projection2D!="latLon" && mag(vectorScale)<VSMALL && maxMagU>VSMALL)
    {
        const scalar meshSize = mag(mesh.bounds().max()-mesh.bounds().min());
        vectorScale = meshSize/(maxMagU*Foam::sqrt(1.*mesh.nCells()));
    }
    else if (projection2D=="latLon"&&mag(vectorScale)<VSMALL && maxMagU>VSMALL)
    {
        vectorScale = 180./(maxMagU*Foam::sqrt(1.*mesh.nCells()));
    }
    Info << "Using vectorScale = " << vectorScale << endl;

    // Set the centre point and valid distance for Orthographic projection
    point p0(0,0,0);
    scalar maxDist2 = 0;
    scalar r = 0;
    if (ortho)
    {
        point pi = mesh.points()[0];
        r = mag(pi);
        p0 = polarPoint(lon0, lat0, r).cartesian(xdir, ydir, axis);
        maxDist2 = 2.5*sqr(r);
    }

    // Find the plot bounding box
    //boundBox bb = setRegion(mesh, projection2D, region, 1.1*r);

    // lon and lat directions if needed
    vectorField rHat = mesh.C();
    vectorField lonHat = axis ^ rHat;
    vectorField latHat = rHat ^ lonHat;
    if (projection2D == "latLon" || projection2D == "latz")
    {
        rHat /= mag(rHat);
        lonHat /= mag(lonHat);
        latHat /= mag(latHat);
    }

    for(label faci = 0; faci < Uf.size(); faci+= vectorFreq)
    {
        point p = mesh.C()[faci];

        // First check that vector for this face should be plot
        bool plotFace = !ortho || magSqr(p - p0) < maxDist2;

        if (plotFace && rMax > rMin)
        {
            const scalar r = mag(p);
            if (r < rMin || r > rMax) plotFace = false;
        }

        if (plotFace && lonMax > lonMin)
        {
            const scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
            if (lon < lonMin || lon > lonMax) plotFace = false;
        }

        if (plotFace && latMax > latMin)
        {
            const scalar lat = radToDeg*Foam::asin(p&axis/mag(p));
            if (lat < latMin || lat > latMax) plotFace = false;
        }

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
            else if (projection2D == "latz")
            {
                scalar r = mag(p);
                const scalar mToDeg = radToDeg/r;
                scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                os << lat << "   " << r << "   "
                   << lat + mToDeg*vectorScale*((Uf[faci]) & (latHat[faci]))<<" "
                   << r + vectorScale*((Uf[faci]) & (rHat[faci])) << '\n';
            }
            else if (projection2D == "lonz")
            {
                scalar r = mag(p);
                const scalar mToDeg = radToDeg/r;
                scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
                os << lon << "   " << r << "   "
                   << lon + mToDeg*vectorScale*((Uf[faci]) & (lonHat[faci]))<<" "
                   << r + vectorScale*((Uf[faci]) & (rHat[faci])) << '\n';
            }
            else if (projection2D == "latLon")
            {
                scalar r = mag(p);
                scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                //scalar sinLat = r > VSMALL ? p.z()/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
                scalar lonEnd = lon + vectorScale*((Uf[faci]) & (lonHat[faci]));
                if (lonEnd >= -180 && lonEnd <= 180 && sign(lon) == sign(lonEnd))
                {
                    os << lon << "   " << lat << "   " << lonEnd << "   "
                       << lat + vectorScale*((Uf[faci]) & (latHat[faci])) << '\n';
                }
            }
            else
            {
                FatalErrorIn("gmtFoam")
                << "projection2D should be one of XY, XZ, YZ, latz, lonz or latLon but "
                    << projection2D << " given"
                    << exit(FatalError);
            }
        }
    }
}

// plot vectors
systemCall = "gmt psxy '" + tmpFile + "' -J -R -h1 -Sv0.3c+s+e"
                  +"+n" + stringScalar(0.3) //100./vectorScale)
                  + " -W" + fieldsToPlot[ifield].colourScale()
                  + " -G" + fieldsToPlot[ifield].colourScale()
                  + " -O -K >> " + epsFile;
systemVerbose(systemCall);

// Add max vector annotation if required
if (writeVectorMax)
{
    systemCall = "echo '" + stringScalar(vectorMaxX) + " "
               + stringScalar(vectorMaxY) + " 12 0 Times-Roman TL "
               + "maximum vector = " + stringScalar(maxMagU,4);
    systemCall += "' | gmt pstext -J -R -h0 -N -O -K >>" + epsFile;
    systemVerbose(systemCall);
}

// remove temp file
systemCall = "rm '" + tmpFile + "' ";
systemVerbose(systemCall);
