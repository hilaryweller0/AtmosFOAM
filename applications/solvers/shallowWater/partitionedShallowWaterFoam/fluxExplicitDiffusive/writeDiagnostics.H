{
    minmax << runTime.timeName() 
           << " " << min(hSum).value() 
           << " " << max(hSum).value()
           << " " << min(sigma[0]).value()
           << " " << max(sigma[0]).value()
           << endl;
           
    total  << runTime.timeName()
           << " " << sum(hSum).value()
           << " " << sum(sigma[0]).value()
           << endl;
}
