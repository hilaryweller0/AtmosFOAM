#include "PolynomialFit.H"
#include "SVD.H"

template<class Polynomial>
Foam::PolynomialFit<Polynomial>::PolynomialFit
(
    const bool linearCorrection,
    const scalar linearLimitFactor,
    const scalar centralWeight, 
    const direction dimensions,
    const label terms
)
:
    linearCorrection_(linearCorrection),
    linearLimitFactor_(linearLimitFactor),
    centralWeight_(centralWeight),
    dim_(dimensions),
    minSize_(terms)
{}

template<class Polynomial>
void Foam::PolynomialFit<Polynomial>::fit
(
    scalarList& coeffsi,
    scalarList& wts,
    const List<point>& C,
    const scalar wLin,
    const point& p0,
    const bool pureUpwind,
    const vector& idir,
    const vector& jdir,
    const vector& kdir,
    const label faceI
)
{
    // Setup the point weights
    wts[0] = centralWeight_;
    if (!pureUpwind)
    {
        wts[1] = centralWeight_;
    }

    // Local coordinate scaling
    scalar scale = 1;

    // Matrix of the polynomial components
    scalarRectangularMatrix B(C.size(), minSize_, scalar(0));

    forAll(C, ip)
    {
        vector d = toLocalCoordinates(p0, C[ip], idir, jdir, kdir);

        if (ip == 0)
        {
            scale = cmptMax(cmptMag((d)));
        }

        // Scale the radius vector
        d /= scale;

        Polynomial::addCoeffs(B[ip], d, wts[ip], dim_);
    }

    // Additional weighting for constant (and linear) terms
    for (label i = 0; i < B.n(); i++)
    {
        B[i][0] *= wts[0];
        B[i][1] *= wts[0];
    }

    // Set the fit
    label stencilSize = C.size();
    coeffsi.setSize(stencilSize);

    bool goodFit = false;
    for (int iIt = 0; iIt < 8 && !goodFit; iIt++)
    {
        SVD svd(B, SMALL);

        scalar maxCoeff = 0;
        label maxCoeffi = 0;

        for (label i=0; i<stencilSize; i++)
        {
            coeffsi[i] = wts[0]*wts[i]*svd.VSinvUt()[0][i];
            if (mag(coeffsi[i]) > maxCoeff)
            {
                maxCoeff = mag(coeffsi[i]);
                maxCoeffi = i;
            }
        }

        if (linearCorrection_)
        {
            goodFit =
                (mag(coeffsi[0] - wLin) < linearLimitFactor_*wLin)
             && (mag(coeffsi[1] - (1 - wLin)) < linearLimitFactor_*(1 - wLin))
             && maxCoeffi <= 1;
        }
        else
        {
            // go through all coeffsi except the first, add up all positive coeffs
            scalar positiveCoeffSum = 0;
            forAll(coeffsi, i)
            {
                if (i > 0 && coeffsi[i] > 0)
                {
                    positiveCoeffSum += coeffsi[i];
                }
            }

            // Upwind: weight on face is 0
            goodFit = (mag(coeffsi[0] - 1.0) < linearLimitFactor_*1.0)
                && coeffsi[0] > positiveCoeffSum;
        }

        if (!goodFit) // (not good fit so increase weight in the centre and
                      //  weight for constant and linear terms)
        {

            wts[0] *= 10;
            if (linearCorrection_)
            {
                wts[1] *= 10;
            }
            else if (!pureUpwind && iIt == 0)
            {
                wts[1] /= centralWeight_;
            }

            for (label j = 0; j < B.m(); j++)
            {
                B[0][j] *= 10;
                if (linearCorrection_) B[1][j] *= 10;
                else if (!pureUpwind && iIt == 0) B[1][j] /= centralWeight_;
            }

            for (label i = 0; i < B.n(); i++)
            {
                B[i][0] *= 10;
                B[i][1] *= 10;
            }
        }
    }

    if (goodFit)
    {
        if (linearCorrection_)
        {
            // Remove the uncorrected linear coefficients
            coeffsi[0] -= wLin;
            coeffsi[1] -= 1 - wLin;
        }
        else
        {
            // Remove the uncorrected upwind coefficients
            coeffsi[0] -= 1.0;
        }
    }
    else
    {
        WarningIn
        (
            "FitData<Polynomial>::calcFit(..)"
        )   << "Could not fit face " << faceI
            << "    Weights = " << coeffsi
            << ", reverting to upwind/linear." << nl
            << "    Linear weights " << wLin << " " << 1 - wLin << endl;

        coeffsi = 0;
        
        if (linearCorrection_)
        {
            coeffsi[0] = 1-wLin;
            coeffsi[1] = -(1-wLin);
        }
    }
}

template<class Polynomial>
vector Foam::PolynomialFit<Polynomial>::toLocalCoordinates(
        const point& origin,
        const point& p,
        const vector& idir,
        const vector& jdir,
        const vector& kdir)
{
    vector d;
    d.x() = (p - origin)&idir;
    d.y() = (p - origin)&jdir;
    #ifndef SPHERICAL_GEOMETRY
    d.z() = (p - origin)&kdir;
    #else
    d.z() = mag(p) - mag(origin);
    #endif
    return d;
}
