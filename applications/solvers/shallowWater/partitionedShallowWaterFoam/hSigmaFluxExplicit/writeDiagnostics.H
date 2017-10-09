{
    minmax << runTime.timeName() 
           << " " << min(hSum).value() 
           << " " << max(hSum).value()
           << " " << min(sigma[0]).value()
           << " " << max(sigma[0]).value()
           << endl;
           
    total  << runTime.timeName()
           << " " << (fvc::domainIntegrate(hSum)/Vtot).value()
           << " " << (fvc::domainIntegrate(sigma[0])/Vtot).value()
           << endl;
}
