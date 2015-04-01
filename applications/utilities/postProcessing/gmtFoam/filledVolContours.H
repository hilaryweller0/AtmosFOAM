// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / vf.name() + ".dat";
{
    OFstream os(tmpFile);
    os << "#x    y    " << vf.name() << endl;

    const fvPatchScalarField& ssf = vf.boundaryField()[patchi];
    Info << vf.name() << " goes from " << min(ssf) << " to " << max(ssf)
         << endl;

    if (colourMin == colourMax)
    {
        colourMin = min(ssf);
        colourMax = max(ssf);
        colourStep = 0.1*(colourMax - colourMin);
    }

    // Set the centre point and valid distance for Orthographic projection
    point p0(0,0,0);
    scalar maxDist2 = 0;
    if (ortho)
    {
        point pi = mesh.points()[plotPatch[0][0]];
        scalar r = mag(pi);
        p0 = polarPoint(lon0, lat0, r).cartesian(xdir, ydir, axis);
        maxDist2 = 2*sqr(r);
    }

    for(label faci = 0; faci < ssf.size(); faci++)
    {
        point p = plotPatch.faceCentres()[faci];

        // First check that this vertex should be plot
        bool plotFace = !ortho || magSqr(p - p0) < maxDist2;

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
            else if (projection2D == "latLon")
            {
                scalar r = mag(p);
                scalar sinLat = r > VSMALL ? (p&axis)/r : 0;
                //scalar sinLat = r > VSMALL ? p.z()/r : 0;
                sinLat = max(min(sinLat, 1), -1);
                scalar lat = radToDeg*Foam::asin(sinLat);
                scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
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
    systemCall += " -D9c/-1.1c/18c/0.5ch -B"
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
                       + stringScalar(2*colourStep) + "  | ps2eps -f -q > " + legNameH;
    systemVerbose(systemCall);

    systemCall = systemCall2 + " -D-0.5c/9c/18c/0.5c -B"
                       + stringScalar(2*colourStep) + "  | ps2eps -f -q > " + legNameV;
    systemVerbose(systemCall);

    chDir("..");
}

// remove temp files
systemCall = "rm colourScale.cpt '" + tmpFile + "' ";
systemVerbose(systemCall);
