/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::cubicUpwindFitPolynomial

Description
    Cubic polynomial for upwind biased interpolation fitting.

    Can be used with the UpwindFit scheme to create a cubic surface
    interpolation scheme

\*---------------------------------------------------------------------------*/

#ifndef cubicUpwindCPCFitPolynomial_H
#define cubicUpwindCPCFitPolynomial_H

#include "vector.H"

namespace Foam
{
    class cubicUpwindCPCFitPolynomial
    {
    public:

        static label nTerms(const direction dim)
        {
            return
            (
                dim == 1 ? 4 :
                dim == 2 ? 12 :
                dim == 3 ? 18 : 0
            );
        }

        static void candidates(List<uint32_t>& candidates, const direction dim)
        {
            if (dim ==1)
            {
                candidates.append(constant | x);
                candidates.append(constant | x | xx);
                candidates.append(constant | x | xx | xxx);
            }
            else if (dim == 2)
            {
                candidates.append(constant | y);
                candidates.append(constant | x);
                candidates.append(constant | y | x);
                candidates.append(constant | y | yy);
                candidates.append(constant | x | xx);
                candidates.append(constant | y | x | xx);
                candidates.append(constant | y | x | xy);
                candidates.append(constant | y | x | yy);
                candidates.append(constant | x | xx | xxx);
                candidates.append(constant | y | x | xx | xxx);
                candidates.append(constant | y | x | xx | xy);
                candidates.append(constant | y | x | xx | yy);
                candidates.append(constant | y | x | xy | yy);
                candidates.append(constant | y | x | xy | yy | xyy);
                candidates.append(constant | y | x | xx | xy | xxx);
                candidates.append(constant | y | x | xx | yy | xxx);
                candidates.append(constant | y | x | xx | xy | xxy);
                candidates.append(constant | y | x | xx | xy | yy);
                candidates.append(constant | y | x | xx | xy | yy | xyy);
                candidates.append(constant | y | x | xx | xy | xxx | xxy);
                candidates.append(constant | y | x | xx | xy | yy | xxx);
                candidates.append(constant | y | x | xx | xy | yy | xxy);
                candidates.append(constant | y | x | xx | xy | yy | xyy | xxx);
                candidates.append(constant | y | x | xx | xy | yy | xyy | xxy);
                candidates.append(constant | y | x | xx | xy | xxx | xxy | xxxy);
                candidates.append(constant | y | x | xx | xy | yy | xxx | xxy);
                candidates.append(constant | y | x | xx | xy | yy | xyy | xxx | xxy);
                candidates.append(constant | y | x | xx | xy | yy | xyy | xxy | xxyy);
                candidates.append(constant | y | x | xx | xy | yy | xxx | xxy | xxxy);
                candidates.append(constant | y | x | xx | xy | yy | xyy | xxx | xxy | xxxy);
                candidates.append(constant | y | x | xx | xy | yy | xyy | xxx | xxy | xxyy);
                candidates.append(constant | y | x | xx | xy | yy | xyy | xxx | xxy | xxxy | xxyy);
                candidates.append(constant | y | x | xx | xy | yy | xyy | xxx | xxy | xxxy | xxyy | xxxyy);
            }
            else if (dim == 3)
            {
                candidates.append(constant | x);
                candidates.append(constant | z);
                candidates.append(constant | y);
                candidates.append(constant | y | yy);
                candidates.append(constant | z | zz);
                candidates.append(constant | x | z);
                candidates.append(constant | x | y);
                candidates.append(constant | x | xx);
                candidates.append(constant | z | y);
                candidates.append(constant | z | y | yz);
                candidates.append(constant | x | y | yy);
                candidates.append(constant | z | y | yy);
                candidates.append(constant | x | y | xy);
                candidates.append(constant | x | z | zz);
                candidates.append(constant | z | y | zz);
                candidates.append(constant | x | z | xz);
                candidates.append(constant | x | z | y);
                candidates.append(constant | x | z | xx);
                candidates.append(constant | x | xx | xxx);
                candidates.append(constant | x | y | xx);
                candidates.append(constant | z | y | yz | yy);
                candidates.append(constant | z | y | yz | zz);
                candidates.append(constant | x | z | y | yz);
                candidates.append(constant | x | y | yy | xy);
                candidates.append(constant | z | y | yy | zz);
                candidates.append(constant | x | z | y | yy);
                candidates.append(constant | x | y | yy | xx);
                candidates.append(constant | x | z | y | xy);
                candidates.append(constant | x | y | xy | xx);
                candidates.append(constant | x | z | zz | xz);
                candidates.append(constant | x | z | y | zz);
                candidates.append(constant | x | z | zz | xx);
                candidates.append(constant | x | z | y | xz);
                candidates.append(constant | x | z | xz | xx);
                candidates.append(constant | x | z | xx | xxx);
                candidates.append(constant | x | z | y | xx);
                candidates.append(constant | x | y | xx | xxx);
                candidates.append(constant | z | y | yz | yy | yyz);
                candidates.append(constant | z | y | yz | yy | zz);
                candidates.append(constant | x | z | y | yz | yy);
                candidates.append(constant | x | z | y | yz | xy);
                candidates.append(constant | z | y | yz | zz | yzz);
                candidates.append(constant | x | z | y | yz | zz);
                candidates.append(constant | x | z | y | yz | xz);
                candidates.append(constant | x | z | y | yz | xx);
                candidates.append(constant | x | y | yy | xy | xyy);
                candidates.append(constant | x | z | y | yy | xy);
                candidates.append(constant | x | y | yy | xy | xx);
                candidates.append(constant | x | z | y | yy | zz);
                candidates.append(constant | x | z | y | yy | xz);
                candidates.append(constant | x | z | y | yy | xx);
                candidates.append(constant | x | y | yy | xx | xxx);
                candidates.append(constant | x | z | y | xy | zz);
                candidates.append(constant | x | y | xy | xx | xxy);
                candidates.append(constant | x | z | y | xy | xz);
                candidates.append(constant | x | z | y | xy | xx);
                candidates.append(constant | x | y | xy | xx | xxx);
                candidates.append(constant | x | z | y | zz | xz);
                candidates.append(constant | x | z | zz | xz | xx);
                candidates.append(constant | x | z | zz | xz | xzz);
                candidates.append(constant | x | z | zz | xx | xxx);
                candidates.append(constant | x | z | y | zz | xx);
                candidates.append(constant | x | z | xz | xx | xxx);
                candidates.append(constant | x | z | xz | xx | xxz);
                candidates.append(constant | x | z | y | xz | xx);
                candidates.append(constant | x | z | y | xx | xxx);
                candidates.append(constant | z | y | yz | yy | zz | yyz);
                candidates.append(constant | x | z | y | yz | yy | yyz);
                candidates.append(constant | x | z | y | yz | yy | xy);
                candidates.append(constant | z | y | yz | yy | zz | yzz);
                candidates.append(constant | x | z | y | yz | yy | zz);
                candidates.append(constant | x | z | y | yz | yy | xz);
                candidates.append(constant | x | z | y | yz | yy | xx);
                candidates.append(constant | x | z | y | yz | xy | zz);
                candidates.append(constant | x | z | y | yz | xy | xz);
                candidates.append(constant | x | z | y | yz | xy | xx);
                candidates.append(constant | x | z | y | yz | zz | yzz);
                candidates.append(constant | x | z | y | yz | zz | xz);
                candidates.append(constant | x | z | y | yz | zz | xx);
                candidates.append(constant | x | z | y | yz | xz | xx);
                candidates.append(constant | x | z | y | yz | xx | xxx);
                candidates.append(constant | x | z | y | yy | xy | zz);
                candidates.append(constant | x | z | y | yy | xy | xyy);
                candidates.append(constant | x | y | yy | xy | xx | xyy);
                candidates.append(constant | x | y | yy | xy | xx | xxy);
                candidates.append(constant | x | z | y | yy | xy | xz);
                candidates.append(constant | x | z | y | yy | xy | xx);
                candidates.append(constant | x | y | yy | xy | xx | xxx);
                candidates.append(constant | x | z | y | yy | zz | xz);
                candidates.append(constant | x | z | y | yy | zz | xx);
                candidates.append(constant | x | z | y | yy | xz | xx);
                candidates.append(constant | x | z | y | yy | xx | xxx);
                candidates.append(constant | x | z | y | xy | zz | xz);
                candidates.append(constant | x | z | y | xy | zz | xx);
                candidates.append(constant | x | z | y | xy | xx | xxy);
                candidates.append(constant | x | y | xy | xx | xxy | xxx);
                candidates.append(constant | x | z | y | xy | xz | xx);
                candidates.append(constant | x | z | y | xy | xx | xxx);
                candidates.append(constant | x | z | zz | xz | xx | xxx);
                candidates.append(constant | x | z | zz | xz | xx | xxz);
                candidates.append(constant | x | z | y | zz | xz | xx);
                candidates.append(constant | x | z | y | zz | xz | xzz);
                candidates.append(constant | x | z | zz | xz | xx | xzz);
                candidates.append(constant | x | z | y | zz | xx | xxx);
                candidates.append(constant | x | z | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | xz | xx | xxx);
                candidates.append(constant | x | z | y | xz | xx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | yyz);
                candidates.append(constant | z | y | yz | yy | zz | yyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | zz | yyz);
                candidates.append(constant | x | z | y | yz | yy | xz | yyz);
                candidates.append(constant | x | z | y | yz | yy | xx | yyz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz);
                candidates.append(constant | x | z | y | yz | yy | xy | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz);
                candidates.append(constant | x | z | y | yz | yy | xy | xx);
                candidates.append(constant | x | z | y | yz | yy | zz | yzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz);
                candidates.append(constant | x | z | y | yz | yy | zz | xx);
                candidates.append(constant | x | z | y | yz | yy | xz | xx);
                candidates.append(constant | x | z | y | yz | yy | xx | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | yzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz);
                candidates.append(constant | x | z | y | yz | xy | zz | xx);
                candidates.append(constant | x | z | y | yz | xy | xz | xyz);
                candidates.append(constant | x | z | y | yz | xy | xx | xxy);
                candidates.append(constant | x | z | y | yz | xy | xz | xx);
                candidates.append(constant | x | z | y | yz | xy | xx | xxx);
                candidates.append(constant | x | z | y | yz | zz | xz | yzz);
                candidates.append(constant | x | z | y | yz | zz | xx | yzz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx);
                candidates.append(constant | x | z | y | yz | zz | xz | xzz);
                candidates.append(constant | x | z | y | yz | zz | xx | xxx);
                candidates.append(constant | x | z | y | yz | xz | xx | xxx);
                candidates.append(constant | x | z | y | yz | xz | xx | xxz);
                candidates.append(constant | x | z | y | yy | xy | zz | xyy);
                candidates.append(constant | x | z | y | yy | xy | zz | xz);
                candidates.append(constant | x | z | y | yy | xy | zz | xx);
                candidates.append(constant | x | y | yy | xy | xx | xyy | xxy);
                candidates.append(constant | x | z | y | yy | xy | xz | xyy);
                candidates.append(constant | x | z | y | yy | xy | xx | xyy);
                candidates.append(constant | x | y | yy | xy | xx | xyy | xxx);
                candidates.append(constant | x | z | y | yy | xy | xx | xxy);
                candidates.append(constant | x | y | yy | xy | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yy | xy | xz | xx);
                candidates.append(constant | x | z | y | yy | xy | xx | xxx);
                candidates.append(constant | x | z | y | yy | zz | xz | xx);
                candidates.append(constant | x | z | y | yy | zz | xz | xzz);
                candidates.append(constant | x | z | y | yy | zz | xx | xxx);
                candidates.append(constant | x | z | y | yy | xz | xx | xxx);
                candidates.append(constant | x | z | y | yy | xz | xx | xxz);
                candidates.append(constant | x | z | y | xy | zz | xx | xxy);
                candidates.append(constant | x | z | y | xy | zz | xz | xx);
                candidates.append(constant | x | z | y | xy | zz | xz | xzz);
                candidates.append(constant | x | z | y | xy | zz | xx | xxx);
                candidates.append(constant | x | z | y | xy | xz | xx | xxy);
                candidates.append(constant | x | z | y | xy | xx | xxy | xxx);
                candidates.append(constant | x | z | y | xy | xz | xx | xxx);
                candidates.append(constant | x | z | y | xy | xz | xx | xxz);
                candidates.append(constant | x | z | zz | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | zz | xz | xx | xxx);
                candidates.append(constant | x | z | zz | xz | xx | xxx | xzz);
                candidates.append(constant | x | z | y | zz | xz | xx | xxz);
                candidates.append(constant | x | z | zz | xz | xx | xxz | xzz);
                candidates.append(constant | x | z | y | zz | xz | xx | xzz);
                candidates.append(constant | x | z | y | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | yyz);
                candidates.append(constant | x | z | y | yz | yy | xy | yyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | yyz);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | yyz);
                candidates.append(constant | x | z | y | yz | yy | zz | yyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | yyz);
                candidates.append(constant | x | z | y | yz | yy | zz | xx | yyz);
                candidates.append(constant | x | z | y | yz | yy | xz | xx | yyz);
                candidates.append(constant | x | z | y | yz | yy | xx | yyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xyz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | xxx);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | yzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xx | yzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xx | xxx);
                candidates.append(constant | x | z | y | yz | yy | xz | xx | xxx);
                candidates.append(constant | x | z | y | yz | yy | xz | xx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xyz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | yzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xx | yzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xx | xxy);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xx | xxx);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xyz);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xxy);
                candidates.append(constant | x | z | y | yz | xy | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xxx);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xxz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | yzz);
                candidates.append(constant | x | z | y | yz | zz | xz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | zz | xx | yzz | xxx);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | xxx);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | xxz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | xzz);
                candidates.append(constant | x | z | y | yz | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xyy);
                candidates.append(constant | x | z | y | yy | xy | zz | xx | xyy);
                candidates.append(constant | x | z | y | yy | xy | zz | xx | xxy);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xx | xxx);
                candidates.append(constant | x | z | y | yy | xy | xx | xyy | xxy);
                candidates.append(constant | x | y | yy | xy | xx | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xyy);
                candidates.append(constant | x | z | y | yy | xy | xx | xyy | xxx);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xxy);
                candidates.append(constant | x | z | y | yy | xy | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xxx);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xxz);
                candidates.append(constant | x | z | y | yy | zz | xz | xx | xxx);
                candidates.append(constant | x | z | y | yy | zz | xz | xx | xxz);
                candidates.append(constant | x | z | y | yy | zz | xz | xx | xzz);
                candidates.append(constant | x | z | y | yy | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxy);
                candidates.append(constant | x | z | y | xy | zz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxx);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xzz);
                candidates.append(constant | x | z | y | xy | xz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | xy | xz | xx | xxy | xxz);
                candidates.append(constant | x | z | y | xy | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | zz | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | zz | xz | xx | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | zz | xz | xx | xxx | xzz);
                candidates.append(constant | x | z | y | zz | xz | xx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | yyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | yyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | yyz | xyz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | yyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | yyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | yyz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | yyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | yyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xx | yyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | yyz | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xx | yyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xz | xx | yyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xz | xx | yyz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xx | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | xxx);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xyz | yzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xyz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xx | yzz | xxy);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xx | yzz | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxy);
                candidates.append(constant | x | z | y | yz | xy | zz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xzz);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xyz | xxy);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xyz | xxx);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xyz | xxz);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xxy | xxz);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | yzz | xxx);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | yzz | xxz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | yzz | xzz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | xxx | xzz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | xxz | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xx | xyy | xxy);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xyy | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xx | xyy | xxx);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxy);
                candidates.append(constant | x | z | y | yy | xy | zz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxx);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xzz);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xyy | xxy);
                candidates.append(constant | x | z | y | yy | xy | xx | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xyy | xxx);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xyy | xxz);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xxy | xxz);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yy | zz | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yy | zz | xz | xx | xxx | xzz);
                candidates.append(constant | x | z | y | yy | zz | xz | xx | xxz | xzz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxy | xxz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxy | xzz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxx | xzz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxz | xzz);
                candidates.append(constant | x | z | y | xy | xz | xx | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | zz | xz | xx | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | yyz | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | yyz | xyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | yyz | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | yyz | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | yyz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | yyz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xx | yyz | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xz | xx | yyz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yzz | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xyz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxy);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxy);
                candidates.append(constant | x | z | y | yz | xy | zz | xx | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxy | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxy | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxx | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xyz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xyz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xyz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxy);
                candidates.append(constant | x | z | y | yy | xy | zz | xx | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxx);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxy | xxz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxy | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxx | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxz | xzz);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xyy | xxy | xxz);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xyy | xxx | xxz);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yy | zz | xz | xx | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xx | yyz | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | yzz | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyz | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyz | xyy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | xyy | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyy | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | xyy | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xyy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyy | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxy);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | xz | xx | xyz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | zz | xz | xx | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxy | xxz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxy | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxx | xxz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxx | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxz | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yy | xy | xz | xx | xyy | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | xy | zz | xz | xx | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyz | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyz | xyy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | xyy | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyy | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | xyy | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xyy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyy | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xyz | xyy | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | xyy | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xyy | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xyy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyy | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | yyz | xyz | xyy | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xx | yyz | xyy | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xyy | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xyy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyy | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | zz | xz | xx | yyz | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | xyz | xyy | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | yzz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yy | xy | zz | xz | xx | xyy | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxy);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | xz | xx | yyz | xyz | xyy | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yzz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | xy | zz | xz | xx | xyz | yzz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxy | xxx);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxy | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxy | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | yzz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | yzz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyy | yzz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxy | xxx | xxz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxy | xxx | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxy | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | yzz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyy | yzz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | xyz | xyy | yzz | xxy | xxx | xxz | xzz);
                candidates.append(constant | x | z | y | yz | yy | xy | zz | xz | xx | yyz | xyz | xyy | yzz | xxy | xxx | xxz | xzz);
            }
        }

        static scalar secondXderivative
        (
            const point& p,
            const uint32_t terms,
            const label stencilIndex,
            const scalarRectangularMatrix& Binv
        )
        {
            scalar dx_sqr = 0.0;
            if (terms & xx) dx_sqr += 2*Binv(rowOf(xx, terms), stencilIndex);
            if (terms & xxx) dx_sqr += 6*p.x()*Binv(rowOf(xxx, terms), stencilIndex);
            if (terms & xxy) dx_sqr += 2*p.y()*Binv(rowOf(xxy, terms), stencilIndex);
            return dx_sqr;
        }

        static void addCoeffs
        (
            scalar* coeffs,
            const vector& d,
            const scalar weight,
            const direction dim
        )
        {
            register label curIdx = 0;

            coeffs[curIdx++] = weight;
            coeffs[curIdx++] = weight*d.x();
            coeffs[curIdx++] = weight*sqr(d.x());
            coeffs[curIdx++] = weight*pow(d.x(),3);
            coeffs[curIdx++] = weight*sqr(d.x())*sqr(d.y());
            coeffs[curIdx++] = weight*pow(d.x(),3)*d.y();
            coeffs[curIdx++] = weight*pow(d.x(),3)*sqr(d.y());

            if (dim >= 2)
            {
                coeffs[curIdx++] = weight*d.y();
                coeffs[curIdx++] = weight*d.x()*d.y();
                coeffs[curIdx++] = weight*sqr(d.y());
                coeffs[curIdx++] = weight*sqr(d.x())*d.y();
                coeffs[curIdx++] = weight*d.x()*sqr(d.y());
            }

            if (dim == 3)
            {
                coeffs[curIdx++] = weight*d.z();
                coeffs[curIdx++] = weight*d.x()*d.z();
                coeffs[curIdx++] = weight*d.y()*d.z();
                coeffs[curIdx++] = weight*sqr(d.z());
                coeffs[curIdx++] = weight*d.x()*d.y()*d.z();
                coeffs[curIdx++] = weight*sqr(d.x())*d.z();
                coeffs[curIdx++] = weight*d.x()*sqr(d.z());
                coeffs[curIdx++] = weight*sqr(d.y())*d.z();
                coeffs[curIdx++] = weight*d.y()*sqr(d.z());
            }
        }

        private:
        static label rowOf(const label term, const uint32_t terms)
        {
            label index = 0;
            for (label i=0; i<=20; i++)
            {
                if ((1 << i) == term && (term & terms)) return index;
                if ((1 << i) && terms) index++;
            }

            return -1;
        }

        static const label constant = 1 << 0;
        static const label x        = 1 << 1;
        static const label xx       = 1 << 2;
        static const label xxx      = 1 << 3;
        static const label y        = 1 << 4;
        static const label xy       = 1 << 5;
        static const label yy       = 1 << 6;
        static const label xxy      = 1 << 7;
        static const label xyy      = 1 << 8;
        static const label xxyy     = 1 << 9;
        static const label xxxy     = 1 << 10;
        static const label xxxyy    = 1 << 11;
        static const label z        = 1 << 12;
        static const label xz       = 1 << 13;
        static const label yz       = 1 << 14;
        static const label zz       = 1 << 15;
        static const label xyz      = 1 << 16;
        static const label xxz      = 1 << 17;
        static const label xzz      = 1 << 18;
        static const label yyz      = 1 << 19;
        static const label yzz      = 1 << 20;
    };
}

#endif
