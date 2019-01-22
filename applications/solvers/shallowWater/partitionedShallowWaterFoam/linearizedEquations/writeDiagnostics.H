{
    minmax << runTime.timeName() 
           << " " << min(h).value() 
           << " " << max(h).value()
           << " " << min(sigma[0]).value()
           << " " << max(sigma[0]).value()
           << endl;
           
    total  << runTime.timeName()
           << " " << (fvc::domainIntegrate(h)/Vtot).value()
           << " " << (fvc::domainIntegrate(sigma[0])/Vtot).value()
           << endl;
}
