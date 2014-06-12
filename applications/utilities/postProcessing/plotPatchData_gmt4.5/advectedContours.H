// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / cf.name() + ".dat";
{
    OFstream os(tmpFile);
    os << "#x    y    " << cf.name() << endl;

    // Set the centre point and valid distance for Orthographic projection
    point p0(0,0,0);
    scalar maxDist2 = 0;
    if (ortho)
    {
        point pi = mesh.points()[plotPatch[0][0]];
        scalar r = mag(pi);
        p0 = polarPoint(lon0, lat0, r).cartesian();
        maxDist2 = 2*sqr(r);
    }

    for(label ci = 0; ci < cf.size(); ci++)
    {
        const contour& c = cf[ci];
        
        // Check that this contour should be plot
        if (c.contourVal() >= colourMin && c.contourVal() <= colourMax)
        {
            os << "> -Z" << c.contourVal() << endl;
        
            forAllConstIter(contour, c, cIter)
//            for(label ip = 0; ip < c.size(); ip++)
            {
                const point& p = cIter().position();
    
                // First check that this point should be plot
                bool plotPoint = !ortho || magSqr(p - p0) < maxDist2;

                if (plotPoint)
                {
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
                        os << lon << "   " << lat << '\n';
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
    }
}

// plot contours
string systemCall = "makecpt -C" + fieldsToPlot[ifield].colourScale();
if (colourStep > 0)
{
    systemCall += " -Z -D -T" + stringScalar(colourMin)/stringScalar(colourMax)
                  /stringScalar(colourStep) + " > colourScale.cpt";
}
else
{
    systemCall += " -I -Z -D -T" +stringScalar(colourMax)/stringScalar(colourMin)
                  /stringScalar(-colourStep) + " > colourScale.cpt";
}
systemVerbose(systemCall);

systemCall = "psxy '" + tmpFile + "' -CcolourScale.cpt"
                           + " -J -R -H1 -m -W+1p -O -K >> " + epsFile;
systemVerbose(systemCall);

systemCall = "psxy '" + tmpFile + "' -J -R -H1 -m -Sp -O -K >> " + epsFile;
systemVerbose(systemCall);

// Plot legend
if (!separateLegend)
{
    systemCall = "psscale -CcolourScale.cpt -S -E -K -O ";
        systemCall += " -D9c/-0.5c/18c/0.5ch -B"
        + stringScalar(2*mag(colourStep)) + " >> " + epsFile;
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

    string systemCall2 = "psscale -C../colourScale.cpt -S -E ";
    systemCall = systemCall2 + " -D9c/-0.5c/18c/0.5ch -B"
               + stringScalar(2*mag(colourStep)) + " > " + legNameH;
    systemVerbose(systemCall);

    systemCall = systemCall2 + " -D-0.5c/9c/18c/0.5c -B"
               + stringScalar(2*mag(colourStep)) + " > " + legNameV;
    systemVerbose(systemCall);

    chDir("..");
}

// remove temp files
systemCall = "rm colourScale.cpt '" + tmpFile + "' ";
systemVerbose(systemCall);
