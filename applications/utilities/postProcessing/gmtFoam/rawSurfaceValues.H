surfaceScalarField ssf(fieldHeader, mesh);

scalar colourMin = fieldsToPlot[ifield].min();
scalar colourMax = fieldsToPlot[ifield].max();
scalar colourStep = fieldsToPlot[ifield].delta();

// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / ssf.name() + ".dat";

bool autoscale = false;

scalar minSsf = 0;
scalar maxSsf = 0;
{
    OFstream os(tmpFile);
    os << "#x    y    " << ssf.name() << endl;

    Info << ssf.name() << " goes from " << min(ssf) << " to " << max(ssf)
         << endl;
    minSsf = min(ssf.internalField());
    maxSsf = max(ssf.internalField());
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
        
        colourStep = (colourMax - colourMin)/22.;
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
        face f = mesh.faces()[faci];

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
        
        // Check that the face centre is within rMin, rMax and lonMin, lonMax
        if (plotFace)
        {
            const point& p = mesh.Cf()[faci];
            if (rMax > rMin)
            {
                const scalar r = mag(p);
                if (r < rMin || r > rMax) plotFace = false;
            }
            if (plotFace && lonMax > lonMin)
            {
                const scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
                if (lon < lonMin || lon > lonMax) plotFace = false;
            }
        }

        if (plotFace)
        {
            // Make face f2 which goes back to the beginning again
            labelList f2(f.size()+1);
            forAll(f,ip) {f2[ip] = f[ip];}
            f2[f.size()] = f[0];

            os << "> -Z" << ssf[faci] << endl;
            forAll (f2, ip)
            {
                point p = mesh.points()[f2[ip]];

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
                    FatalErrorIn("gmtFoam")
                   << "projection2D should be one of XY, XZ, YZ or latLon but "
                        << projection2D << " given" << exit(FatalError);
                }
            }
        }
    }
}

// colourscale for plot
systemCall = "makecpt -C" + fieldsToPlot[ifield].colourScale()
                + " -Z -D -T" + stringScalar(colourMin)/stringScalar(colourMax)
                           /stringScalar(colourStep) + " > colourScale.cpt";
systemVerbose(systemCall);

// main colour plot
systemCall = "psxy '" + tmpFile + "' -CcolourScale.cpt"
                           + " -J -R -h1 -W2pt -O -K >> " + epsFile;
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
    if (!autoscale)
    {
        systemCall2 = systemCall2 + " -B" + stringScalar(2*mag(colourStep));
    }
    systemCall = systemCall2 + " -D9c/-0.5c/18c/0.5ch  | ps2eps -f -q > " + legNameH;
    systemVerbose(systemCall);

    systemCall = systemCall2 + " -D-0.5c/9c/18c/0.5c -B"
                       + stringScalar(2*colourStep) + "  | ps2eps -f -q > " + legNameV;
    systemVerbose(systemCall);

    chDir("..");
}

// Add max and min text if required
if (annotate == "minMax")
{
    systemCall = "echo '" + stringScalar(annotateX) + " "
               + stringScalar(annotateY) + " 12 0 Times-Roman TL "
               + "min = " + stringScalar(minSsf,5)
               + "    max = " + stringScalar(maxSsf,5);
//               + "    L@-1@- = " + stringScalar(l1,4)
//               + "    L@-2@- = " + stringScalar(l2,4);
    systemCall += "' | pstext -J -R -h0 -N -O -K >>" + epsFile;
    systemVerbose(systemCall);
}


// delete temporary files
systemCall = "rm colourScale.cpt '" + tmpFile + "' ";
systemVerbose(systemCall);
