// The FOAM Project // File: testSperical.C
/*
-------------------------------------------------------------------------------
 =========         | Application
 \\      /         |
  \\    /          | Name:   testSperical
   \\  /           | Family: spherical
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION
    test the spherical geometry classes

AUTHOR
    Hilary Spencer.

-------------------------------------------------------------------------------
*/ 

#include "polarPoint.H"
#include "polarVector.H"
#include "lCartPolarVector.H"
#include "physicalConstants.H"

using namespace Foam;
using namespace physicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    static const scalar earthRadius = 6371e3;
    
    polarPoint pp(0, 45*pi/180, earthRadius);
    point pc = convertToCart(pp);
    polarPoint pp2 = convertToPolar(pc);

    lCartPolarVector lcv(40,0,0);
    
    polarVector pv = convertToPolar(lcv, pc);
    polarVector pv2 = convertToPolar(lcv, pp);
        
    vector v = convertToCart(lcv, pc);
    vector v2 = convertToCart(lcv, pp);
    vector v3 = convertToCart(pv, pc);
    vector v4 = convertToCart(pv, pp);
    
//    polarVector pv3 = convertToPolar(v, pc);
//    polarVector pv4 = convertToPolar(v, pp);
    
//    lCartPolarVector lcv1 = convertToLCart(v, pc);
//    lCartPolarVector lcv2 = convertToLCart(v, pp);
    lCartPolarVector lcv3 = convertToLCart(pv, pc);
    lCartPolarVector lcv4 = convertToLCart(pv, pp);
    
    
    Info << "pp = " << pp << " pc = " << pc << " pp2 = " << pp2 << endl;
    
    Info << "lcv = " << lcv 
         << " lcv3 = " << lcv3 << " lcv4 = " << lcv4 << endl;
         
    Info << "v = " << v << " v2 = " << v2 << " v3 = " << v3 << " v4 = " << v4 
         << endl;
    Info << "pv = " << pv << " pv2 = " << pv2 << endl;

    Info << "End\n" << endl;
    
    return 0;
}


// ************************************************************************* //
