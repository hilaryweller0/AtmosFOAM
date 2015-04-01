if (annotate != "" && annotate != "minMax")
{
    if (annotate == "day")
    {
        const scalar t = runTime.value();
        label nDays = floor(t/86400.);
        label nHours = round((t - nDays*86400)/3600.);
        systemCall = "echo '" + stringScalar(annotateX) + " "
                   + stringScalar(annotateY) + " 12 0 4 TL "
                   + stringScalar(nDays) + " days";
        if (nHours) systemCall += ", " + stringScalar(nHours) + " hours";

    }
    else
    {
        systemCall = "echo '" + stringScalar(annotateX) + " "
                   + stringScalar(annotateY) + " 12 0 Times-Roman TL " + annotate;
    }
    
    systemCall += "' | pstext -J -R -h0 -N -K -O >>" + epsFile;
    systemVerbose(systemCall);
}

if (xlabel != "")
{
    systemCall = "echo '" + stringScalar(xlabelLoc.first()) + " "
               + stringScalar(xlabelLoc.second()) + " 12 0 Times-Roman TL " + xlabel;
    
    systemCall += "' | pstext -J -R -h0 -N -K -O >>" + epsFile;
    systemVerbose(systemCall);
}

if (ylabel != "")
{
    systemCall = "echo '" + stringScalar(ylabelLoc.first()) + " "
               + stringScalar(ylabelLoc.second()) + " 12 90 Times-Roman TL " + ylabel;
    
    systemCall += "' | pstext -J -R -h0 -N -K -O >>" + epsFile;
    systemVerbose(systemCall);
}

