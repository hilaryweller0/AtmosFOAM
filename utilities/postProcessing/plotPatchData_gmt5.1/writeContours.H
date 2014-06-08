// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / vf.name() + ".dat";

{
    OFstream os(tmpFile);
    os << "#x    y    " << vf.name() << endl;
    
    pointField centres(mesh.nCells() + mesh.nFaces() - mesh.nInternalFaces());
    scalarField vals(centres.size());
    label nPoints = mesh.nCells();
    label ip = mesh.nCells();
    forAll(vf, cellI)
    {
        centres[cellI] = mesh.C()[cellI];
        vals[cellI] = vf[cellI];
    }
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];
        if (vf.boundaryField()[patchI].size() != 0)
        {
            nPoints += patch.size();
            forAll(patch, facei)
            {
                centres[ip] = patch.faceCentres()[facei];
                vals[ip] = vf.boundaryField()[patchI][facei];
                ip++;
            }
        }
    }
    vals.setSize(nPoints);
    centres.setSize(nPoints);

    Info << vf.name() << " goes from " << min(vals) << " to "
         << max(vals) << endl;
    if (colourMin == colourMax)
    {
        colourMin = min(vals);
        colourMax = max(vals);
        colourStep = 0.1*(colourMax - colourMin);
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


    // Plot the values
    for(label ip = 0; ip < vals.size(); ip++)
    {
        point p = centres[ip];

        // First check that this point should be plot
        bool plotPoint = !ortho || magSqr(p - p0) < maxDist2;

        if (plotPoint && rMax > rMin)
        {
            const scalar r = mag(p);
            if (r < rMin || r > rMax) plotPoint = false;
        }

        if (plotPoint && lonMax > lonMin)
        {
            const scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
            if (lon < lonMin || lon > lonMax) plotPoint = false;
        }

        if (plotPoint && latMax > latMin)
        {
            const scalar lat = radToDeg*Foam::asin(p&axis/mag(p));
            if (lat < latMin || lat > latMax) plotPoint = false;
        }

        if (plotPoint)
        {
            if (projection2D == "XY")
            {
                os << p.x() << "   " << p.y() << "   " << vals[ip] << '\n';
            }
            else if (projection2D == "YX")
            {
                os << p.y() << "   " << p.x() << "   " << vals[ip] << '\n';
            }
            else if (projection2D == "XZ")
            {
                os << p.x() << "   " << p.z() << "   " << vals[ip] << '\n';
            }
            else if (projection2D == "ZX")
            {
                os << p.z() << "   " << p.x() << "   " << vals[ip] << '\n';
            }
            else if (projection2D == "YZ")
            {
                os << p.y() << "   " << p.z() << "   " << vals[ip] << '\n';
            }
            else if (projection2D == "ZY")
            {
                os << p.z() << "   " << p.y() << "   " << vals[ip] << '\n';
            }
            else if (projection2D == "latz")
            {
                scalar r = mag(p);
                scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                os << lat << "   " << r << "   " << vals[ip] << '\n';
            }
            else if (projection2D == "lonz")
            {
                scalar r = mag(p);
                scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
                os << lon << "   " << r << "   " << vals[ip] << '\n';
            }
            else if (projection2D == "latLon")
            {
                scalar r = mag(p);
                scalar sinLat = r > VSMALL ? p.z()/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                scalar lon = radToDeg*Foam::atan2(p.y(), p.x());
                os << lon << "   " << lat << "   " << vals[ip] << '\n';
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
           + " -A- -J -R -h1 -W" + contourType + " -D -O -K >> " + epsFile;
systemVerbose(systemCall);

// remove temp files
//systemCall = "rm contours.cpt '" + tmpFile + "' ";
//systemVerbose(systemCall);

