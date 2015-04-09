surfaceScalarField ssf(fieldHeader, mesh);

// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / ssf.name() + ".dat";
{
    OFstream os(tmpFile);
    os << "#x    y    " << ssf.name() << endl;

    Info << ssf.name() << " goes from " << min(ssf) << " to " << max(ssf)
         << endl;
    if (mag(colourStep) <= VSMALL)
    {
        colourMin = min(ssf).value();
        colourMax = max(ssf).value();
        colourStep = 0.1*(colourMax - colourMin);
    }
    if (colourMax-colourMin < 5e-4)
    {
        colourMin -= 1;
        colourMax += 1;
        colourStep = 1;
    }

    // Set the centre point and valid distance for Orthographic projection
    point p0(0,0,0);
    scalar maxDist2 = 0;
    if (ortho)
    {
        point pi = mesh.points()[0];
        scalar r = mag(pi);
        p0 = polarPoint(lon0, lat0, r).cartesian();
        maxDist2 = 2*sqr(r);
    }

    for(label faci = 0; faci < mesh.nInternalFaces(); faci++)
    {
        point p = mesh.faceCentres()[faci];

        // First check that this point should be plot
        bool plotPoint = !ortho || magSqr(p - p0) < maxDist2;

        if (plotPoint)
        {
            if (projection2D == "XY")
            {
                os << p.x() << "   " << p.y() << "   " << ssf[faci] << '\n';
            }
            else if (projection2D == "YX")
            {
                os << p.y() << "   " << p.x() << "   " << ssf[faci] << '\n';
            }
            else if (projection2D == "XZ")
            {
                os << p.x() << "   " << p.z() << "   " << ssf[faci] << '\n';
            }
            else if (projection2D == "ZX")
            {
                os << p.z() << "   " << p.x() << "   " << ssf[faci] << '\n';
            }
            else if (projection2D == "YZ")
            {
                os << p.y() << "   " << p.z() << "   " << ssf[faci] << '\n';
            }
            else if (projection2D == "ZY")
            {
                os << p.z() << "   " << p.y() << "   " << ssf[faci] << '\n';
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
                os << lon << "   " << lat << "   " << ssf[faci] << '\n';
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

// plot contours
systemCall = "makecpt -N -T"
         + stringScalar(colourMin)/stringScalar(colourMax)
           /stringScalar(colourStep) + " > contours.cpt";
systemVerbose(systemCall);

string contourType = fieldsToPlot[ifield].colourScale();
if(fieldsToPlot[ifield].plotType() == FieldToPlot::DASHED_CONTOURS)
{
    contourType = "-";
}

systemCall = "pscontour '" + tmpFile + "' -Ccontours.cpt"
           + " -A- -J -R -h1 -W" + contourType + " -O -K >> " + epsFile;
systemVerbose(systemCall);

// remove temp files
systemCall = "rm contours.cpt '" + tmpFile + "' ";
systemVerbose(systemCall);
