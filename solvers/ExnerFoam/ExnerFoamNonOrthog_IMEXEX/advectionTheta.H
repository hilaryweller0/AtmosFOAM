// Sub step for advecting theta using and RK3 scheme
{
    dimensionedScalar dtsub = dt/nAdvectionSubSteps;

    Info << "Theta advection sub steps" << flush;
    for(int subStep = 0; subStep < nAdvectionSubSteps; subStep++)
    {
        Info << " " << subStep << flush;
        divUt1 = fvc::div(U,theta);
        theta = (theta.oldTime()*rho.oldTime() - dtsub*divUt1)/rho;
        
        divUt2 = fvc::div(U,theta);
        theta = (theta.oldTime()*rho.oldTime() - dtsub*(0.75*divUt1 - 0.25*divUt2))/rho;
        
        divUt3 = fvc::div(U,theta);
        theta = (theta.oldTime()*rho.oldTime() - dtsub*(1/12.*(divUt1 + divUt2) - 2/3.*divUt3))/rho;
    }
    
    Info << endl;
}
