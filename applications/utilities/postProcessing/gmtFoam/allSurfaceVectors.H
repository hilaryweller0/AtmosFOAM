// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.name() / U.name() + ".dat";
scalar maxMagU = 0;
{
    OFstream os(tmpFile);
    os << "#x    y    " << U.name() << endl;

    const DimensionedField<vector,surfaceMesh>& Uf = U.internalField();
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
    if (ortho)
    {
        point pi = 1.01*mesh.Cf()[0];
        scalar r = mag(pi);
        p0 = polarPoint(lon0, lat0, r).cartesian(xdir, ydir, axis);
        maxDist2 = 2*sqr(r);
    }
    // Find the plot bounding box
    boundBox bb = setRegion(mesh, projection2D, region, 1.01*mag(mesh.Cf()[0]));

    // lon and lat directions if needed
    vectorField rHat = mesh.Cf();
    vectorField lonHat = axis ^ rHat;
    vectorField latHat = rHat ^ lonHat;
    if (projection2D == "latLon")
    {
        //Info << "lonHat = " << lonHat << " latHat = " << latHat << nl;    
        rHat /= max(mag(rHat), VSMALL);
        lonHat /= max(mag(lonHat), VSMALL);
        latHat /= max(mag(latHat), VSMALL);
    }

    for(label faci = 0; faci < Uf.size(); faci+= vectorFreq)
    {
        point p = mesh.Cf()[faci];

        // First check that vector for this face should be plot
        bool plotFace = (!ortho || magSqr(p - p0) < maxDist2)
                      && bb.contains(p);
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
                os << p.x() - 0.5*vectorScale*Uf[faci][0] << "   "
                   << p.y() - 0.5*vectorScale*Uf[faci][1] << "   "
                   << p.x() + 0.5*vectorScale*Uf[faci][0] << "   "
                   << p.y() + 0.5*vectorScale*Uf[faci][1] << '\n';
            }
            else if (projection2D == "YX")
            {
                os << p.y() - 0.5*vectorScale*Uf[faci][1] << "   " 
                   << p.x() - 0.5*vectorScale*Uf[faci][0] << "   "
                   << p.y() + 0.5*vectorScale*Uf[faci][1] << "   "
                   << p.x() + 0.5*vectorScale*Uf[faci][0] << '\n';
            }
            else if (projection2D == "XZ")
            {
                os << p.x() - 0.5*vectorScale*Uf[faci][0] << "   "
                   << p.z() - 0.5*vectorScale*Uf[faci][2] << "   "
                   << p.x() + 0.5*vectorScale*Uf[faci][0] << "   "
                   << p.z() + 0.5*vectorScale*Uf[faci][2] << '\n';
            }
            else if (projection2D == "ZX")
            {
                os << p.z() - 0.5*vectorScale*Uf[faci][2] << "   "
                   << p.x() - 0.5*vectorScale*Uf[faci][0] << "   "
                   << p.z() + 0.5*vectorScale*Uf[faci][2] << "   "
                   << p.x() + 0.5*vectorScale*Uf[faci][0] << '\n';
            }
            else if (projection2D == "YZ")
            {
                os << p.y() - 0.5*vectorScale*Uf[faci][1] << "   "
                   << p.z() - 0.5*vectorScale*Uf[faci][2] << "   "
                   << p.y() + 0.5*vectorScale*Uf[faci][1] << "   "
                   << p.z() + 0.5*vectorScale*Uf[faci][2] << '\n';
            }
            else if (projection2D == "ZY")
            {
                os << p.z() - 0.5*vectorScale*Uf[faci][2] << "   "
                   << p.y() - 0.5*vectorScale*Uf[faci][1] << "   "
                   << p.z() + 0.5*vectorScale*Uf[faci][2] << "   "
                   << p.y() + 0.5*vectorScale*Uf[faci][1] << '\n';
            }
            else if (projection2D == "latz")
            {
                scalar r = mag(p);
                const scalar mToDeg = radToDeg/r;
                scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                os << lat - 0.5*mToDeg*vectorScale*((Uf[faci]) & (latHat[faci]))<<"   "
                   << r - 0.5*vectorScale*((Uf[faci]) & (rHat[faci])) << "   "
                   << lat + 0.5*mToDeg*vectorScale*((Uf[faci]) & (latHat[faci]))<<" "
                   << r + 0.5*vectorScale*((Uf[faci]) & (rHat[faci])) << '\n';
            }
            else if (projection2D == "lonz")
            {
                scalar r = mag(p);
                const scalar mToDeg = radToDeg/r;
                scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
                os << lon - 0.5*mToDeg*vectorScale*((Uf[faci]) & (lonHat[faci]))<<"   "
                   << r - 0.5*vectorScale*((Uf[faci]) & (rHat[faci])) << "   "
                   << lon + 0.5*mToDeg*vectorScale*((Uf[faci]) & (lonHat[faci]))<<" "
                   << r + 0.5*vectorScale*((Uf[faci]) & (rHat[faci])) << '\n';
            }
            else if (projection2D == "latLon")
            {
                scalar r = mag(p);
                scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                //scalar sinLat = r > VSMALL ? p.z()/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat)
                           - 0.5*vectorScale*((Uf[faci]) & (latHat[faci]));
                scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir))
                           - 0.5*vectorScale*((Uf[faci]) & (lonHat[faci]));
                scalar lon2 = lon + vectorScale*((Uf[faci]) & (lonHat[faci]));
                scalar lat2 = lat + vectorScale*((Uf[faci]) & (latHat[faci]));
                if
                (
                    sign(lon2 - lonMin) == sign(lon - lonMin)
                 && lon2 < 180 && lon2 > -180 && lon < 180 && lon > -180
                )
                {
                    os << lon << "   " << lat << "   " << lon2 << "   " << lat2
                       << '\n';
                }
            }
            else
            {
                FatalErrorIn("gmtFoam")
                << "projection2D should be one of XY, XZ, YZ, latz or latLon but "
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
