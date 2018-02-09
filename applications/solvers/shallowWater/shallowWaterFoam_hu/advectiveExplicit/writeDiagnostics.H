{
    minmax << runTime.timeName() 
           << " " << min(h).value() 
           << " " << max(h).value()
           << endl;
           
    total  << runTime.timeName()
           << " " << (fvc::domainIntegrate(h)/Vtot).value()
           << endl;
}
