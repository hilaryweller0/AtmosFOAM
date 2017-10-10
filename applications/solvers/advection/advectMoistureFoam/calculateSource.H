forAll(Sv, celli)
{
    //Condense vapour into liquid.
    if (rv[celli] > rvs[celli])      
    { 
        Sl[celli] = 0;
        Sv[celli] = 1./(timeScale*dt.value());
    }
    //Evaporate liquid into vapour.
    else if (rv[celli] < rvs[celli])
    {
        //If there is surplus liquid, evaporate to equilibrium.
        if ( (rvs[celli]-rv[celli]) <= rl[celli] )    
        {
            Sl[celli] = 0;
            Sv[celli] = 1./(timeScale*dt.value());
        }
        //Evaporate all the liquid available.
        else                                          
        {
            Sl[celli] = 1./(timeScale*dt.value());
            Sv[celli] = 0;
        }
    }
    //If rv=rvs, ensure Sl and Sv are set to 0 to avoid artefacts.
    else
    { 
        Sl[celli] = 0;
        Sv[celli] = 0;
    }
}
