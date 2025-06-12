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
        
        colourStep = (colourMax - colourMin)/22.;
        if (colourMax-colourMin < 5e-6)
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

    for(label celli = 0; celli < ssf.size(); celli++)
    {
        // Cell centre
        const point& C = mesh.C()[celli];
        const polarPoint Cpp = (projection2D == "latLon") ? 
                               convertToPolar(C, lonMax) : polarPoint();
    
        // First check that all vertices of this cell should be plot
        bool plotCell = true;

        for(label faci = 0; faci < mesh.cells()[celli].size()&&plotCell;faci++)
        {
            face f = mesh.faces()[mesh.cells()[celli][faci]];

            // Check that all vertices of this face should be plot
            if (ortho)
            {
                for(label ip = 0; ip < f.size() && plotCell; ip++)
                {
                    point p = mesh.points()[f[ip]];
                    if (magSqr(p - p0) >= maxDist2) plotCell = false;
                    if (!bb.contains(p)) plotCell = false;
                }
            }
            else if(projection2D == "cylinder")
            {
                scalar thetaPrev = 0;
                for(label ip = 0; ip < f.size() && plotCell; ip++)
                {
                    point p = mesh.points()[f[ip]];
                    scalar maxTheta = 360;
                    if (region != "g") maxTheta = bb.max()[0];
                    cylinderPoint cp = radToDeg*convertToCylinder(p,maxTheta);
                    if (ip == 0) thetaPrev = cp.theta();
                    if
                    (
                        cp.theta() < bb.min().x()
                     || cp.theta() > bb.max().x()
                     || cp.z() < bb.min().z()
                     || cp.z() > bb.max().z()
                     || mag(cp.theta() - thetaPrev) > 180
                    )
                    {
                        plotCell = false;
                    }
                    thetaPrev = cp.theta();
                }
            }
        }
        
        // Check that the cell centre is within rMin, rMax and lonMin, lonMax
        if (plotCell && rMax > rMin)
        {
            const scalar r = mag(C);
            if (r < rMin || r > rMax) plotCell = false;
        }
        if (plotCell && lonMax > lonMin)
        {
            const point& p = mesh.C()[celli];
            const scalar lon = radToDeg*Foam::atan2((p&ydir), (p&xdir));
            if (lon < lonMin || lon > lonMax) plotCell = false;
        }
        
        if (plotCell)
        {
            if (!plotZeroValuePolys && mag(ssf[celli]) < 1e-14) continue;
            for(label faci = 0; faci < mesh.cells()[celli].size(); faci++)
            {
                bool plotFace = true;
                // For lat-lon projection, check that this is not a vertical face
                if (plotCell && projection2D == "latLon")
                {
                    label faceI = mesh.cells()[celli][faci];
                    if
                    (
                        (mesh.faceCentres()[faceI] & mesh.faceAreas()[faceI])
                < 0.2*mag(mesh.faceCentres()[faceI])*mag(mesh.faceAreas()[faceI])
                    )
                    {
                        plotFace = false;
                    }
                }
                
                if (plotFace)
                {
                os << "> -Z" << ssf[celli] << endl;

                face f = mesh.faces()[mesh.cells()[celli][faci]];

                forAll (f, ip)
                {
                    point p = mesh.points()[f[ip]];

                    if (projection2D == "XY") os << p.x()<<"   "<< p.y()<<'\n';
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
                    else if (projection2D == "cylinder")
                    {
                        scalar maxTheta = 360;
                        if (region != "g") maxTheta = bb.max()[0];
                        cylinderPoint cp = radToDeg*convertToCylinder(p,maxTheta);
                        os << cp.theta() << " " << cp.z() <<nl;
                    }
                    else if (projection2D == "latLon")
                    {
/*                        scalar r = mag(p);
                        scalar sinLat = r > VSMALL ? p.z()/r : 0;
                        sinLat = max(min(sinLat, 1), -1);
                        scalar lat = radToDeg*Foam::asin(sinLat);
                        scalar cosLat = Foam::sqrt(1 - sqr(sinLat));
                        scalar cosLon = mag(r*cosLat) > VSMALL ?
                                p.x()/(r*cosLat):0;
                        cosLon = max(min(cosLon, 1), -1);
                        scalar lon = radToDeg*sign(p.y())*Foam::acos(cosLon);
*/
                        polarPoint ppp = convertToPolar(p, lonMax);
                        scalar lat = polarPoint::radToDeg*ppp.lat();
                        scalar lon = polarPoint::radToDeg*ppp.lon();
                        
                        if
                        (
                            lonMin != lonMax &&
                            mag(polarPoint::radToDeg*Cpp.lon() - lon) < 60
                        )
                        {
                            os << lon << "   " << lat << '\n';
                        }
                    }
                    else
                    {
                        FatalErrorIn("gmtFoam")
                   << "projection2D should be one of XY, XZ, YZ, cylinder or latLon but "
                       << projection2D << " given" << exit(FatalError);
                    }
                }
                }
            }
        }
    }
}

// colourscale for plot
systemCall = "gmt makecpt -C" + fieldsToPlot[ifield].colourScale();
if (colourStep > 0)
{
    systemCall += " -Z -D -T" + stringScalar(colourMin,6)/stringScalar(colourMax,6)
                  /stringScalar(colourStep,6) + " > colourScale.cpt";
}
else
{
    systemCall += " -I -Z -D -T" +stringScalar(colourMax,4)/stringScalar(colourMin,4)
                  /stringScalar(-colourStep,3) + " > colourScale.cpt";
}
systemVerbose(systemCall);

// main colour plot
systemCall = "gmt psxy '" + tmpFile + "' -CcolourScale.cpt"
                           + " -J -R -h1 -L -O -K >> " + epsFile;
systemVerbose(systemCall);

// Plot legend
if (!separateLegend)
{
    systemCall = "gmt psscale -CcolourScale.cpt -S -K -O ";
    //if (!autoscale)
    {
        systemCall = systemCall + " -B" + stringScalar(2*mag(colourStep),3);
    }
    systemCall = systemCall + " -D9c/-0.5c/18c/0.5ch >> " + epsFile;
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
    systemCall += "' | gmt pstext -J -R -h0 -N -O -K >>" + epsFile;
    systemVerbose(systemCall);
}

// delete temporary files
systemCall = "rm colourScale.cpt '" + tmpFile + "' ";
systemVerbose(systemCall);
