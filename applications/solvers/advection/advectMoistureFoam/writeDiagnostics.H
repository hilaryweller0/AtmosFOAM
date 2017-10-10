{
    //Minima and maxima of fields.
    minmax_prognostic << runTime.timeName() 
                      << " " << min(rl).value() 
                      << " " << max(rl).value()
                      << " " << min(rv).value()
                      << " " << max(rv).value()
                      << " " << min(rv+rl).value()
                      << " " << max(rv+rl).value()
                      << endl;
                      
    minmax_diagnostic << runTime.timeName() 
                      << " " << min(rl_diag).value() 
                      << " " << max(rl_diag).value()
                      << " " << min(rv_diag).value()
                      << " " << max(rv_diag).value()
                      << " " << min(rt).value()
                      << " " << max(rt).value()
                      << endl;
                      
    minmax_analytic   << runTime.timeName() 
                      << " " << min(rl_analytic).value() 
                      << " " << max(rl_analytic).value()
                      << " " << min(rv_analytic).value()
                      << " " << max(rv_analytic).value()
                      << " " << min(rt_analytic).value()
                      << " " << max(rt_analytic).value()
                      << endl;
    
    //Conservation properties of fields.           
    rl_total = (fvc::domainIntegrate(rl)/Vtot).value();
    rv_total = (fvc::domainIntegrate(rv)/Vtot).value();
    rt_total = (fvc::domainIntegrate(rl+rv)/Vtot).value();
    total_prognostic  << runTime.timeName()
                      << " " << rl_total
                      << " " << rv_total
                      << " " << rt_total
                      << endl;
                      
    rl_total = (fvc::domainIntegrate(rl_diag)/Vtot).value();
    rv_total = (fvc::domainIntegrate(rv_diag)/Vtot).value();
    rt_total = (fvc::domainIntegrate(rt)/Vtot).value();
    total_diagnostic  << runTime.timeName()
                      << " " << rl_total
                      << " " << rv_total
                      << " " << rt_total
                      << endl;
                      
    rl_total = (fvc::domainIntegrate(rl_analytic)/Vtot).value();
    rv_total = (fvc::domainIntegrate(rv_analytic)/Vtot).value();
    rt_total = (fvc::domainIntegrate(rt_analytic)/Vtot).value();
    total_analytic  << runTime.timeName()
                      << " " << rl_total
                      << " " << rv_total
                      << " " << rt_total
                      << endl;
                      
    //l2 normalised error relative to analytic solution.
    rt2_total = fvc::domainIntegrate( sqr(rt_analytic) ).value();
    error_prognostic  << runTime.timeName()
                      << " " << Foam::sqrt(fvc::domainIntegrate( sqr(rl - rl_analytic) ).value()/rt2_total)
                      << " " << Foam::sqrt(fvc::domainIntegrate( sqr(rv - rv_analytic) ).value()/rt2_total)
                      << endl;
                      
    error_diagnostic  << runTime.timeName()
                      << " " << Foam::sqrt(fvc::domainIntegrate( sqr(rl_diag - rl_analytic) ).value()/rt2_total)
                      << " " << Foam::sqrt(fvc::domainIntegrate( sqr(rv_diag - rv_analytic) ).value()/rt2_total)
                      << endl;
}
