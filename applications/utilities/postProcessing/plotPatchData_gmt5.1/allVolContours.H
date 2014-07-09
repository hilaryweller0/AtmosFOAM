// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / vf.name() + ".dat";

// File for a boundary line
fileName tmpBFile = args.rootPath() / args.caseName()
                   / runTime.timeName() / vf.name() +boundaryPatchName+".dat";
{
    OFstream os(tmpFile);
    os << "#x    y    " << vf.name() << endl;

    const scalarField& sf = vf.internalField();
    Info << vf.name() << " goes from " << min(sf) << " to " << max(sf)
        << endl;
    if (colourMin == colourMax)
    {
        colourMin = min(sf);
        colourMax = max(sf);
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

    // Plot a boundary patch if requested
    const label bpatchi = boundaryPatchName == "" ? -1 :
                          mesh.boundaryMesh().findPatchID(boundaryPatchName);

    const polyPatch& plotBPatch = mesh.boundaryMesh()[max(bpatchi, 0)];

    if (bpatchi != -1)
    {
        OFstream osb(tmpBFile);

        const fvPatchScalarField& ssf = vf.boundaryField()[bpatchi];

        for
        (
            label faci = plotBPatch.start();
            faci < plotBPatch.start() + plotBPatch.size();
            faci++
        )
        {
            point p = mesh.faceCentres()[faci];
            label i = faci - plotBPatch.start();

            // First check that this point should be plot
            bool plotPoint = !ortho || magSqr(p - p0) < maxDist2;

            if (plotPoint)
            {
                //point pOut = p + 1e-3*plotBPatch.faceNormals()[faci];
                //point pOut =p+1e-6*plotBPatch[faci].normal(plotBPatch.points());
                point pOut = p - 0.01*(mesh.C()[mesh.faceOwner()[faci]] - p);
                point pOut2 = p - 0.02*(mesh.C()[mesh.faceOwner()[faci]] - p);

                if (projection2D == "XY")
                {
                    os << p.x() << "   " << p.y() << "   " << ssf[i] << '\n';
                    //os << pOut.x()<<"   "<<pOut.y()<< "   NaN" << '\n';
                    //os << pOut2.x()<<"   "<<pOut2.y()<< "   NaN" << '\n';
                }
                else if (projection2D == "YX")
                {
                    os << p.y() << "   " << p.x() << "   " << ssf[i] << '\n';
                    //os << pOut.y()<<"   "<<pOut.x()<< "   NaN" << '\n';
                    //os << pOut2.y()<<"   "<<pOut2.x()<< "   NaN" << '\n';
                }
                else if (projection2D == "XZ")
                {
                    os << p.x() << "   " << p.z() << "   " << ssf[i] << '\n';
                    //os << pOut.x()<<"   "<<pOut.z()<< "   NaN" << '\n';
                    //os << pOut2.x()<<"   "<<pOut2.z()<< "   NaN" << '\n';
                }
                else if (projection2D == "ZX")
                {
                    os << p.z() << "   " << p.x() << "   " << ssf[i] << '\n';
                    //os << pOut.z()<<"   "<<pOut.x()<< "   NaN" << '\n';
                }
                else if (projection2D == "YZ")
                {
                    os << p.y() << "   " << p.z() << "   " << ssf[i] << '\n';
                    //os << pOut.y()<<"   "<<pOut.z()<< "   NaN" << '\n';
                    //os << pOut2.y()<<"   "<<pOut2.z()<< "   NaN" << '\n';
                }
                else if (projection2D == "ZY")
                {
                    os << p.z() << "   " << p.y() << "   " << ssf[i] << '\n';
                    os << pOut.z()<<"   "<<pOut.y()<< "   NaN" << '\n';
                    os << pOut2.z()<<"   "<<pOut2.y()<< "   NaN" << '\n';
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
                    os << lon << "   " << lat << "   " << ssf[i] << '\n';

                    scalar sinLatOut =  r > VSMALL ? pOut.z()/r : 0;
                    sinLatOut = max(min(sinLatOut, 1), -1);
                    scalar latOut = radToDeg*Foam::asin(sinLatOut);
                    scalar cosLatOut = Foam::sqrt(1 - sqr(sinLatOut));
                    scalar cosLonOut = mag(r*cosLatOut) > VSMALL ?
                                                   pOut.x()/(r*cosLatOut):0;
                    cosLonOut = max(min(cosLonOut, 1), -1);
                    scalar lonOut =radToDeg*sign(pOut.y())*Foam::acos(cosLonOut);
                    os << lonOut << "   " << latOut << "   inf\n";
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

    // Plot the internal values
    for(label faci = 0; faci < sf.size(); faci++)
    {
        point p = mesh.C()[faci];

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
                os << p.x() << "   " << p.y() << "   " << sf[faci] << '\n';
            }
            else if (projection2D == "YX")
            {
                os << p.y() << "   " << p.x() << "   " << sf[faci] << '\n';
            }
            else if (projection2D == "XZ")
            {
                os << p.x() << "   " << p.z() << "   " << sf[faci] << '\n';
            }
            else if (projection2D == "ZX")
            {
                os << p.z() << "   " << p.x() << "   " << sf[faci] << '\n';
            }
            else if (projection2D == "YZ")
            {
                os << p.y() << "   " << p.z() << "   " << sf[faci] << '\n';
            }
            else if (projection2D == "ZY")
            {
                os << p.z() << "   " << p.y() << "   " << sf[faci] << '\n';
            }
            else if (projection2D == "latz")
            {
                scalar r = mag(p);
                scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                os << lat << "   " << r << "   " << sf[faci] << '\n';
            }
            else if (projection2D == "lonz")
            {
                scalar r = mag(p);
                scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
                os << lon << "   " << r << "   " << sf[faci] << '\n';
            }
            else if (projection2D == "latLon")
            {
                scalar r = mag(p);
                scalar sinLat = r > VSMALL ? p.z()/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                scalar lon = radToDeg*Foam::atan2(p.y(), p.x());
                os << lon << "   " << lat << "   " << sf[faci] << '\n';
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
           + " -A- -J -R -h1 -W" + contourType + " -O -K >> " + epsFile;
systemVerbose(systemCall);

