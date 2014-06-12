if (annotate != "" && annotate != "minMax")
{
    string systemCall = "mv .gmt_bb_info .gmt_bb_infoSave";
    systemVerbose(systemCall);
    systemCall = "mv .gmtcommands4 .gmtcommands4Save";
    systemVerbose(systemCall);

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
    
    systemCall += "' | pstext -JX15c -R0/1/0/1 -H0 -N -K -O >>" + epsFile;
    systemVerbose(systemCall);

    systemCall = "mv .gmt_bb_infoSave .gmt_bb_info";
    systemVerbose(systemCall);
    systemCall = "mv .gmtcommands4Save .gmtcommands4";
    systemVerbose(systemCall);
}

if (xlabel != "")
{
    string systemCall = "mv .gmt_bb_info .gmt_bb_infoSave";
    systemVerbose(systemCall);
    systemCall = "mv .gmtcommands4 .gmtcommands4Save";
    systemVerbose(systemCall);

    systemCall = "echo '" + stringScalar(xlabelLoc.first()) + " "
               + stringScalar(xlabelLoc.second()) + " 12 0 Times-Roman TL " + xlabel;
    
    systemCall += "' | pstext -JX15c -R0/1/0/1 -H0 -N -K -O >>" + epsFile;
    systemVerbose(systemCall);

    systemCall = "mv .gmt_bb_infoSave .gmt_bb_info";
    systemVerbose(systemCall);
    systemCall = "mv .gmtcommands4Save .gmtcommands4";
    systemVerbose(systemCall);
}

if (ylabel != "")
{
    string systemCall = "mv .gmt_bb_info .gmt_bb_infoSave";
    systemVerbose(systemCall);
    systemCall = "mv .gmtcommands4 .gmtcommands4Save";
    systemVerbose(systemCall);

    systemCall = "echo '" + stringScalar(ylabelLoc.first()) + " "
               + stringScalar(ylabelLoc.second()) + " 12 90 Times-Roman TL " + ylabel;
    
    systemCall += "' | pstext -JX15c -R0/1/0/1 -H0 -N -K -O >>" + epsFile;
    systemVerbose(systemCall);

    systemCall = "mv .gmt_bb_infoSave .gmt_bb_info";
    systemVerbose(systemCall);
    systemCall = "mv .gmtcommands4Save .gmtcommands4";
    systemVerbose(systemCall);
}

