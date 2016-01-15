fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.timeName() / vf.name() + ".boundary.dat";

scalar minSsf = VGREAT;
scalar maxSsf = -VGREAT;
{
    OFstream os(tmpFile);
    os << "#x    y    " << vf.name() << " boundary" << endl;


    forAll(vf.boundaryField(), patchI)
    {
        fvPatchField<scalar>& patch = vf.boundaryField()[patchI];
        minSsf = min(minSsf, min(patch));
        maxSsf = max(maxSsf, max(patch));
    }

    Info << vf.name() << " boundary goes from " << minSsf << " to " << maxSsf <<endl;

    if (colourMin == colourMax)
    {
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
        if (colourMax-colourMin < 5e-4)
        {
            colourMin -= 1;
            colourMax += 1;
            colourStep = 1;
        }
    }

    // Set the centre point and valid distance for Orthographic projection
    point p0(0,0,0);
    scalar r = 0;
    if (ortho)
    {
        point pi = mesh.points()[0];
        r = mag(pi);
        p0 = polarPoint(lon0, lat0, r).cartesian(xdir, ydir, axis);
    }

    // Find the plot bounding box
    boundBox bb = setRegion(mesh, projection2D, region, 1.1*r);

    // iterate over boundary faces and plot their values
    forAll(vf.boundaryField(), patchI)
    {
        fvPatchField<scalar>& patch = vf.boundaryField()[patchI];
        forAll(patch, faceI)
        {
            label patchOffset = mesh.boundary()[patchI].start();
            face f = mesh.faces()[patchOffset + faceI];
            labelList f2(f.size()+1);
            forAll(f,ip) {f2[ip] = f[ip];}
            f2[f.size()] = f[0];

            os << "> -Z" << patch[faceI] << endl;
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
                else if (projection2D == "cylinder")
                {
                    scalar maxTheta = 360;
                    if (region != "g") maxTheta = bb.max()[0];
                    cylinderPoint cp = radToDeg*convertToCylinder(p, maxTheta);
                    os << cp.theta() << " " << cp.z() << nl;
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
                   << "projection2D should be one of XY, XZ, YZ, cylinder or latLon but "
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
                           + " -J -R -h1 -W2p -O -K >> " + epsFile;
systemVerbose(systemCall);

// delete temporary files
systemCall = "rm colourScale.cpt '" + tmpFile + "' ";
systemVerbose(systemCall);
