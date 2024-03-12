pointMesh2 pMesh(mesh);
pointScalarField psf(fieldHeader, pMesh);

scalar colourMin = fieldsToPlot[ifield].min();
scalar colourMax = fieldsToPlot[ifield].max();
scalar colourStep = fieldsToPlot[ifield].delta();

// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / psf.name() + ".dat";
{
    OFstream os(tmpFile);
    os << "#x    y    " << psf.name() << endl;

    const DimensionedField<scalar,pointMesh>& ppsf
         = psf.boundaryField()[patchi].internalField();
    Info << psf.name() << " goes from " << min(ppsf).value() << " to "
         << max(ppsf).value() << endl;
    if (colourMin == colourMax)
    {
        colourMin = min(ppsf).value();
        colourMax = max(ppsf).value();
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
    scalar r = 0;
    if (ortho)
    {
        point pi = mesh.points()[plotPatch[0][0]];
        r = mag(pi);
        p0 = polarPoint(lon0, lat0, r).cartesian();
        maxDist2 = 2*sqr(r);
    }

    // Find the plot bounding box
    boundBox bb = setRegion(mesh, projection2D, region, 1.1*r);

    for(label ip = 0; ip < ppsf.size(); ip++)
    {
        point p = plotPatch.points()[ip];

        // First check that this vertex should be plot
        bool plotFace = !ortho || magSqr(p - p0) < maxDist2;

        if (plotFace)
        {
            if (projection2D == "XY")
            {
                os << p.x() << "   " << p.y() << "   " << ppsf[ip] << '\n';
            }
            else if (projection2D == "YX")
            {
                os << p.y() << "   " << p.x() << "   " << ppsf[ip] << '\n';
            }
            else if (projection2D == "XZ")
            {
                os << p.x() << "   " << p.z() << "   " << ppsf[ip] << '\n';
            }
            else if (projection2D == "ZX")
            {
                os << p.z() << "   " << p.x() << "   " << ppsf[ip] << '\n';
            }
            else if (projection2D == "YZ")
            {
                os << p.y() << "   " << p.z() << "   " << ppsf[ip] << '\n';
            }
            else if (projection2D == "ZY")
            {
                os << p.z() << "   " << p.y() << "   " << ppsf[ip] << '\n';
            }
            else if (projection2D == "cylinder")
            {
                scalar maxTheta = 360;
                if (region != "g") maxTheta = bb.max()[0];
                cylinderPoint cp = radToDeg*convertToCylinder(p, maxTheta);
                os << cp.theta() << " " << cp.z() << " "<<ppsf[ip]<<nl;
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
                os << lon << "   " << lat << "   " << ppsf[ip] << '\n';
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

// colourscale for plot
systemCall = "gmt makecpt -D -C" + fieldsToPlot[ifield].colourScale()
                   + " -D -T" + stringScalar(colourMin)/stringScalar(colourMax)
                       /stringScalar(colourStep) + " > colourScale.cpt";
systemVerbose(systemCall);

// plot filled contours
systemCall = "gmt pscontour '" + tmpFile + "' -CcolourScale.cpt"
                       + " -A- -J -R -h1 -I -O -K >> " + epsFile;
systemVerbose(systemCall);

// Plot legend
if (!separateLegend)
{
    systemCall = "gmt psscale -CcolourScale.cpt -S -K -O ";
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

    systemCall2 = "gmt psscale -C../colourScale.cpt -S ";
    systemCall = systemCall2 + " -D9c/-0.5c/18c/0.5ch -B"
                       + stringScalar(2*colourStep) + "  | ps2eps -f -q > " + legNameH;
    systemVerbose(systemCall);

    systemCall = systemCall2 + " -D-0.5c/9c/18c/0.5c -B"
                       + stringScalar(2*colourStep) + "  | ps2eps -f -q > " + legNameV;
    systemVerbose(systemCall);

    chDir("..");
}

// remove temp files
systemCall = "rm colourScale '" + tmpFile + "' ";
systemVerbose(systemCall);
