// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / vf.name() + ".dat";
        
scalar minVf = 0;
scalar maxVf = 0;
{
    OFstream os(tmpFile);
    os << "#x    y    " << vf.name() << endl;

    const scalarField& ssf = vf.internalField();
    minVf = min(vf).value();
    maxVf = max(vf).value();
    Info << vf.name() << " goes from " << minVf
         << " to " << maxVf << endl;

    if (colourMin == colourMax)
    {
        colourMin = minVf;
        colourMax = maxVf;
        colourStep = 0.1*(colourMax - colourMin);
    }
    if (colourMax-colourMin < 5e-6)
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
        p0 = polarPoint(lon0, lat0, r).cartesian(xdir, ydir, axis);
        maxDist2 = 2*sqr(r);
    }

    // Plot a boundary patch if requested
    const label bpatchi = boundaryPatchName == "" ? -1 :
                          mesh.boundaryMesh().findPatchID(boundaryPatchName);

    const polyPatch& plotBPatch = mesh.boundaryMesh()[max(bpatchi, 0)];

    if (bpatchi != -1)
    {
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
                    os << pOut.x()<<"   "<<pOut.y()<< "   NaN" << '\n';
                    os << pOut2.x()<<"   "<<pOut2.y()<< "   NaN" << '\n';
                }
                else if (projection2D == "YX")
                {
                    os << p.y() << "   " << p.x() << "   " << ssf[i] << '\n';
                    os << pOut.y()<<"   "<<pOut.x()<< "   NaN" << '\n';
                    os << pOut2.y()<<"   "<<pOut2.x()<< "   NaN" << '\n';
                }
                else if (projection2D == "XZ")
                {
                    os << p.x() << "   " << p.z() << "   " << ssf[i] << '\n';
                    os << pOut.x()<<"   "<<pOut.z()<< "   NaN" << '\n';
                    os << pOut2.x()<<"   "<<pOut2.z()<< "   NaN" << '\n';
                }
                else if (projection2D == "ZX")
                {
                    os << p.z() << "   " << p.x() << "   " << ssf[i] << '\n';
                    os << pOut.z()<<"   "<<pOut.x()<< "   NaN" << '\n';
                }
                else if (projection2D == "YZ")
                {
                    os << p.y() << "   " << p.z() << "   " << ssf[i] << '\n';
                    os << pOut.y()<<"   "<<pOut.z()<< "   NaN" << '\n';
                    os << pOut2.y()<<"   "<<pOut2.z()<< "   NaN" << '\n';
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
                    FatalErrorIn("gmtFoam")
                    << "projection2D should be one of XY, XZ, YZ or latLon but "
                        << projection2D << " given"
                        << exit(FatalError);
                }
            }
        }
    }

    for(label faci = 0; faci < ssf.size(); faci++)
    {
        point p = mesh.C()[faci];

        // First check that this vertex should be plot
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
            else if (projection2D == "latz")
            {
                scalar r = mag(p);
                scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                os << lat << "   " << r << "   " << ssf[faci] << '\n';
            }
            else if (projection2D == "lonz")
            {
                scalar r = mag(p);
                scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
                os << lon << "   " << r << "   " << ssf[faci] << '\n';
            }
            else if (projection2D == "latLon")
            {
                scalar r = mag(p);
                scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                //scalar sinLat = r > VSMALL ? p.z()/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
                //scalar lon = radToDeg*Foam::atan2(p.y(), p.x());
                os << lon << "   " << lat << "   " << ssf[faci] << '\n';
            }
            else
            {
                FatalErrorIn("gmtFoam")
                << "projection2D should be one of XY, XZ, YZ, lonz, latz or latLon but "
                    << projection2D << " given"
                    << exit(FatalError);
            }
        }
    }
}

// colourscale for plot
systemCall = "makecpt -D -C" + fieldsToPlot[ifield].colourScale()
                  + " -T" + stringScalar(colourMin)/stringScalar(colourMax)
                            /stringScalar(colourStep) + " > colourScale.cpt";
systemVerbose(systemCall);

// plot filled contours
systemCall = "pscontour '" + tmpFile + "' -CcolourScale.cpt"
                       + " -A- -J -R -h1 -I -O -K >> " + epsFile;
systemVerbose(systemCall);

// Plot legend
if (!separateLegend)
{
    systemCall = "psscale -CcolourScale.cpt -S -E -K -O ";
    systemCall += " -D9c/-0.5c/18c/0.5ch -B"
               + stringScalar(2*colourStep) + " >> " + epsFile;
    systemVerbose(systemCall);
}
else
{
    mkDir("legends");
    chDir("legends");

    fileName legNameH = epsFileName;
    if (epsFileName != fieldHeader.name()) legNameH += "_" + fieldHeader.name();
    fileName legNameV = legNameH;
    legNameH += ".eps";
    legNameV += "_v.eps";
    Info << "Creating legends/" << legNameH << " and legends/" << legNameV
        << endl;

    systemCall2 = "psscale -C../colourScale.cpt -S -E ";
    systemCall = systemCall2 + " -D9c/-0.5c/18c/0.5ch -B"
                       + stringScalar(4*colourStep) + " | ps2eps -f -q > " + legNameH;
    systemVerbose(systemCall);

    systemCall = systemCall2 + " -D-0.5c/9c/18c/0.5c -B"
                       + stringScalar(4*colourStep) + " | ps2eps -f -q > " + legNameV;
    systemVerbose(systemCall);

    chDir("..");
}

// Add max and min text if required
if (annotate == "minMax")
{
    systemCall = "echo '" + stringScalar(annotateX) + " "
               + stringScalar(annotateY) + " 12 0 Times-Roman TL "
               + "min = " + stringScalar(minVf,5)
               + "    max = " + stringScalar(maxVf,5);
//               + "    L@-1@- = " + stringScalar(l1,4)
//               + "    L@-2@- = " + stringScalar(l2,4);
    systemCall += "' | pstext -J -R -h0 -N -O -K >>" + epsFile;
    systemVerbose(systemCall);
}

// remove temp files
systemCall = "rm colourScale.cpt '" + tmpFile + "' ";
systemVerbose(systemCall);
