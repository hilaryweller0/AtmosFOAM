// write data to temporary output file
fileName tmpFile = args.rootPath() / args.caseName()
        / runTime.name() / U.name() + ".dat";
{
    OFstream os(tmpFile);
    os << "#x    y    " << U.name() << endl;

    Info << U.name() << " goes from " << min(mag(U)) << " to " << max(mag(U))
        << endl;

    for(label ip = 0; ip < U.size(); ip++)
    {
        if (projection2D == "XY")
        {
            os << U[ip][0] << "   " << U[ip][1] << '\n';
        }
        else if (projection2D == "YX")
        {
            os << U[ip][1] << "   " << U[ip][0] << '\n';
        }
        else if (projection2D == "XZ")
        {
            os << U[ip][0] << "   " << U[ip][2] << '\n';
        }
        else if (projection2D == "ZX")
        {
            os << U[ip][2] << "   " << U[ip][0] << '\n';
        }
        else if (projection2D == "YZ")
        {
            os << U[ip][1] << "   " << U[ip][2] << '\n';
        }
        else if (projection2D == "ZY")
        {
            os << U[ip][2] << "   " << U[ip][1] << '\n';
        }
        else if (projection2D == "latLon")
        {
            scalar r = mag(U[ip]);
            scalar sinLat = r > VSMALL ? U[ip].z()/r : 0;
            sinLat = max(min(sinLat, 1), -1);
            scalar lat = radToDeg*Foam::asin(sinLat);
            scalar cosLat = Foam::sqrt(1 - sqr(sinLat));
            scalar cosLon = mag(r*cosLat) > VSMALL ? U[ip].x()/(r*cosLat):0;
            cosLon = max(min(cosLon, 1), -1);
            scalar lon = radToDeg*sign(U[ip].y())*Foam::acos(cosLon);

            os << lon << "   " << lat << '\n';
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

// plot vectors
systemCall = "gmt psxy '" + tmpFile + "' -J -R -h1 -Sx0.1i"
                  + " -W -O -K >> " + epsFile;
systemVerbose(systemCall);

// remove temp file
systemCall = "rm '" + tmpFile + "' ";
systemVerbose(systemCall);
