// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / vf.name() + ".dat";
{
    OFstream os(tmpFile);
    os << "#x    y    " << vf.name() << endl;

    const fvPatchScalarField& ssf = vf.boundaryField()[patchi];
    Info << vf.name() << " on patch " << plotPatch.name() << " goes from "
        << min(ssf) << " to " << max(ssf) << endl;
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
        p0 = polarPoint(lon0, lat0, r).cartesian();
        maxDist2 = 2*sqr(r);
    }

    for(label faci = 0; faci < ssf.size(); faci++)
    {
        face f = plotPatch[faci];

        // First check that all vertices of this face should be plot
        bool plotFace = true;
        if (ortho)
        {
            for(label ip = 0; ip < f.size() && plotFace; ip++)
            {
                point p = mesh.points()[f[ip]];
                if (magSqr(p - p0) >= maxDist2) plotFace = false;
            }
        }

        if (plotFace)
        {
            os << "> -Z" << ssf[faci] << endl;
            forAll (f, ip)
            {
                point p = mesh.points()[f[ip]];

                if (projection2D == "XY") os << p.x() << "   "<< p.y() << '\n';
                else if (projection2D == "YX")
                    os << p.y() << "   "<< p.x() << '\n';
                else if (projection2D == "XZ")
                    os << p.x() << "   "<< p.z() << '\n';
                else if (projection2D == "ZX")
                    os << p.z() << "   "<< p.x() << '\n';
                else if (projection2D == "YZ")
                    os << p.y() << "   "<< p.z() << '\n';
                else if (projection2D == "ZY")
                    os << p.z() << "   "<< p.y() << '\n';
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
                        << projection2D << " given" << exit(FatalError);
                }
            }
        }
    }
}

// colourscale for plot
systemCall = "makecpt -C" + fieldsToPlot[ifield].colourScale();
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
/*if (handEditColourScale != "false" && handEditColourScale != "none")
{
    if (handEditColourScale == "true")
    {
        systemCall = "vi colourScale.cpt";
        systemVerbose(systemCall);
    }
    else
    {
        systemCall = handEditColourScale + " colourScale.cpt";
        systemVerbose(systemCall);
    }
}
*/
// main colour plot
systemCall = "psxy '" + tmpFile + "' -CcolourScale.cpt"
                           + " -J -R -h1 -L -O -K >> " + epsFile;
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

    systemCall2 = "psscale -C../colourScale.cpt -S -E ";
    systemCall = systemCall2 + " -D9c/-0.5c/18c/0.5ch -B"
               + stringScalar(2*mag(colourStep)) + "  | ps2eps -f -q > " + legNameH;
    systemVerbose(systemCall);

    systemCall = systemCall2 + " -D-0.5c/9c/18c/0.5c -B"
               + stringScalar(2*mag(colourStep)) + "  | ps2eps -f -q > " + legNameV;
    systemVerbose(systemCall);

    chDir("..");
}

// delete temporary files
systemCall = "rm colourScale.cpt '" + tmpFile + "' ";
systemVerbose(systemCall);
