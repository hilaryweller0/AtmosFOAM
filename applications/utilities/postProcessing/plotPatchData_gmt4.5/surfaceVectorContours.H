scalar colourMin = fieldsToPlot[ifield].min();
scalar colourMax = fieldsToPlot[ifield].max();
scalar colourStep = fieldsToPlot[ifield].delta();
const int vectorDir = fieldsToPlot[ifield].vectorDir();

// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / U.name() + ".dat";
{
    OFstream os(tmpFile);
    os << "#x    y    " << U.name() << endl;

    scalarField Ufd = vectorDir >= 0 && vectorDir <= 2 ?
                      (U.internalField()).component(vectorDir) :
                      mag(U.internalField());

    Info << U.name() << '[' << vectorDir << "] goes from " << min(Ufd)
         << " to " << max(Ufd) << endl;
    if (mag(colourStep) <= VSMALL)
    {
        colourMin = min(Ufd);
        colourMax = max(Ufd);
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

    for(label faci = 0; faci < mesh.nInternalFaces(); faci++)
    {
        point p = mesh.faceCentres()[faci];

        // First check that this point should be plot
        bool plotPoint = (!ortho || magSqr(p - p0) < maxDist2);

        if (plotPoint)
        {
            if (projection2D == "XY")
            {
                os << p.x() << "   " << p.y() << "   " << Ufd[faci] << '\n';
            }
            else if (projection2D == "YX")
            {
                os << p.y() << "   " << p.x() << "   " << Ufd[faci] << '\n';
            }
            else if (projection2D == "XZ")
            {
                os << p.x() << "   " << p.z() << "   " << Ufd[faci] << '\n';
            }
            else if (projection2D == "ZX")
            {
                os << p.z() << "   " << p.x() << "   " << Ufd[faci] << '\n';
            }
            else if (projection2D == "YZ")
            {
                os << p.y() << "   " << p.z() << "   " << Ufd[faci] << '\n';
            }
            else if (projection2D == "ZY")
            {
                os << p.z() << "   " << p.y() << "   " << Ufd[faci] << '\n';
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
                os << lon << "   " << lat << "   " << Ufd[faci] << '\n';
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
string systemCall = "makecpt -N -T"
         + stringScalar(colourMin)/stringScalar(colourMax)
           /stringScalar(colourStep) + " > contours.cpt";
systemVerbose(systemCall);

string contourType = fieldsToPlot[ifield].colourScale();
if(fieldsToPlot[ifield].plotType() == FieldToPlot::DASHED_CONTOURS)
{
    contourType = "-";
}

systemCall = "pscontour '" + tmpFile + "' -Ccontours.cpt"
           //+ " -A+kblack -J -R -H1 -W"
           + " -A- -J -R -H1 -W"
           + contourType + " -O -K >> " + epsFile;
systemVerbose(systemCall);

// remove temp files
systemCall = "rm contours.cpt '" + tmpFile + "' ";
systemVerbose(systemCall);
