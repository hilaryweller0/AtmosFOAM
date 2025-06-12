// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.name() / vf.name() + ".dat";
bool autoscale = false;

scalar minSsf = 0;
scalar maxSsf = 0;
{
    OFstream os(tmpFile);
    os << "#x    y    " << vf.name() << endl;

    const DimensionedField<scalar,volMesh>& ssf = vf.internalField();
    Info << vf.name() << " goes from " << min(ssf).value() << " to "
         << max(ssf).value() << endl;

    systemCall = "echo '%% " + vf.name() + " min_max "
               + stringScalar(min(ssf).value()) + " "
               + stringScalar(max(ssf).value()) + "' >> " + epsFile;
    systemVerbose(systemCall);

    minSsf = min(ssf).value();
    maxSsf = max(ssf).value();
    if (colourMin == colourMax)
    {
        autoscale = true;
        
        if (minSsf < 0 && maxSsf > 0)
        {
            colourMin = min(minSsf, -maxSsf);
            colourMax = -colourMin;
        }
        else
        {
            colourMin = minSsf;
            colourMax = maxSsf;
        }
        colourStep = (colourMax - colourMin)/21.;
        if (colourMax-colourMin < 1e-9)
        {
            colourMin -= 1;
            colourMax += 1;
            colourStep = 1;
        }
    }

    // Set the centre point and valid distance for Orthographic projection
    point p0(0,0,0);
    scalar maxDist2 = 0;
    point pi = mesh.C()[0];

    scalar r = mag(pi);
    if (ortho)
    {
        p0 = polarPoint(lon0, lat0, r).cartesian();
        maxDist2 = 2.5*sqr(r);
    }
    // Find the plot bounding box
    boundBox bb = setRegion(mesh, projection2D, region, 1.1*r);

    // Plot a boundary patch if requested
    const label bpatchi = boundaryPatchName == "" ? -1 :
                          mesh.boundaryMesh().findIndex(boundaryPatchName);

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
                else if (projection2D == "cylinder")
                {
                    scalar maxTheta = 360;
                    if (region != "g") maxTheta = bb.max()[0];

                    cylinderPoint cp = convertToCylinder(p,maxTheta);
                    cylinderPoint cpOut = convertToCylinder(pOut,maxTheta);
                    cylinderPoint cpOut2 = convertToCylinder(pOut2,maxTheta);
                    os << radToDeg*cp.theta() << " " << cp.z() << " "
                       << ssf[faci] << nl;
                    os << radToDeg*cpOut.theta() << " " << cpOut.z() << " NaN\n";
                    os << radToDeg*cpOut2.theta() << " " << cpOut2.z()<<" NaN\n";
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
                    << "projection2D should be one of XY, XZ, YZ, cylinder or latLon but "
                        << projection2D << " given"
                        << exit(FatalError);
                }
            }
        }
    }

    for(label celli = 0; celli < ssf.size(); celli++)
    {
        const point& p = mesh.C()[celli];

        // First check that this cell should be plot
        bool plotCell = true;

        if (ortho)
        {
            if (magSqr(p - p0) >= maxDist2) plotCell = false;
            if (!bb.contains(p)) plotCell = false;
        }
        else if(projection2D == "cylinder")
        {
            scalar maxTheta = 360;
            if (region != "g") maxTheta = bb.max()[0];
            cylinderPoint cp = radToDeg*convertToCylinder(p,maxTheta);
            if
            (
                cp.theta() < bb.min().x()
             || cp.theta() > bb.max().x()
             || cp.z() < bb.min().z()
             || cp.z() > bb.max().z()
            )
            {
                plotCell = false;
            }
        }
        
        // Check that the cell centre is within rMin, rMax and lonMin, lonMax
        if (plotCell && rMax > rMin)
        {
            const scalar r = mag(p);
            if (r < rMin || r > rMax) plotCell = false;
        }
        if (plotCell && lonMax > lonMin)
        {
            const scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
            if (lon < lonMin || lon > lonMax) plotCell = false;
        }
        
        if (plotCell)
        {
            if (projection2D == "XY")
            {
                os << p.x() << "   " << p.y() << "   " << ssf[celli] << '\n';
            }
            else if (projection2D == "YX")
            {
                os << p.y() << "   " << p.x() << "   " << ssf[celli] << '\n';
            }
            else if (projection2D == "XZ")
            {
                os << p.x() << "   " << p.z() << "   " << ssf[celli] << '\n';
            }
            else if (projection2D == "ZX")
            {
                os << p.z() << "   " << p.x() << "   " << ssf[celli] << '\n';
            }
            else if (projection2D == "YZ")
            {
                os << p.y() << "   " << p.z() << "   " << ssf[celli] << '\n';
            }
            else if (projection2D == "ZY")
            {
                os << p.z() << "   " << p.y() << "   " << ssf[celli] << '\n';
            }
            else if (projection2D == "latz")
            {
                scalar r = mag(p);
                scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                os << lat << "   " << r << "   " << ssf[celli] << '\n';
            }
            else if (projection2D == "lonz")
            {
                scalar r = mag(p);
                scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
                os << lon << "   " << r << "   " << ssf[celli] << '\n';
            }
            else if (projection2D == "cylinder")
            {
                scalar maxTheta = 360;
                if (region != "g") maxTheta = bb.max()[0];
                cylinderPoint cp = radToDeg*convertToCylinder(p, maxTheta);
                os << cp.theta() << " " << cp.z() << " "<<ssf[celli]<<nl;
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
                os << lon << "   " << lat << "   " << ssf[celli] << '\n';
            }
            else
            {
                FatalErrorIn("gmtFoam")
                << "projection2D should be one of XY, XZ, YZ, lonz, latz, cylinder or latLon but "
                    << projection2D << " given"
                    << exit(FatalError);
            }
        }
    }
}

// colourscale for plot
systemCall = "gmt makecpt -C" + fieldsToPlot[ifield].colourScale();
if (colourStep > 0)
{
    systemCall += " -D -T" + stringScalar(colourMin,6)/stringScalar(colourMax,6)
                  /stringScalar(colourStep,6) + " > colourScale.cpt";
}
else
{
    systemCall += " -I -D -T" +stringScalar(colourMax,4)/stringScalar(colourMin,4)
                  /stringScalar(-colourStep,3) + " > colourScale.cpt";
}
systemVerbose(systemCall);

// plot filled contours
systemCall = "gmt pscontour '" + tmpFile + "' -CcolourScale.cpt"
                       + " -Sp -A- -J -R -h1 -I -O -K >> " + epsFile;
systemVerbose(systemCall);

// Plot legend
if (!separateLegend)
{
    systemCall = "gmt psscale -CcolourScale.cpt -S -K -O ";
    systemCall += " -D9c/-0.5c/18c/0.5ch -B"
               + stringScalar(2*mag(colourStep),3) + " >> " + epsFile;
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

    systemCall2 = "gmt psscale -C../colourScale.cpt -S ";
    if (!autoscale)
    {
        systemCall2 = systemCall2 + " -B" + stringScalar(2*mag(colourStep));
    }
    systemCall = systemCall2 + " -D9c/-0.5c/18c/0.5ch  | ps2eps -f -q > " + legNameH;
    systemVerbose(systemCall);

    systemCall = systemCall2 + " -D-0.5c/9c/18c/0.5c "
                + " | ps2eps -f -q > " + legNameV;
    systemVerbose(systemCall);

    chDir("..");
}

// Add max and min text if required
if (annotate == "minMax")
{
    systemCall = "echo '" + stringScalar(annotateX) + " "
               + stringScalar(annotateY) + " 12 0 Times-Roman TL "
               + "min = " + stringScalar(minSsf,6)
               + "    max = " + stringScalar(maxSsf,6);
//               + "    L@-1@- = " + stringScalar(l1,4)
//               + "    L@-2@- = " + stringScalar(l2,4);
    systemCall += "' | gmt pstext -J -R -h0 -N -O -K >>" + epsFile;
    systemVerbose(systemCall);
}

// remove temp files
systemCall = "rm colourScale.cpt '" + tmpFile + "' ";
systemVerbose(systemCall);
