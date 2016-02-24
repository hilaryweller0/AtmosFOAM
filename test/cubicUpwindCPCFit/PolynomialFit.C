#include "PolynomialFit.H"

Test::PolynomialFit::PolynomialFit()
{
    coefficients_.setSize(12);
    coefficients_[0] = 0.7;
}

const Foam::scalarList Test::PolynomialFit::coefficients()
{
    return coefficients_;
}
