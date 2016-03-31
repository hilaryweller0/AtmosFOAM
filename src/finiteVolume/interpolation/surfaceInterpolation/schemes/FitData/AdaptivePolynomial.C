#include "AdaptivePolynomial.H"
#include "SVD.H"
#include "ListOps.H"

template<class Polynomial>
Foam::AdaptivePolynomial<Polynomial>::AdaptivePolynomial(
        const localStencil& stencil,
        const direction dimensions,
        const scalar tolerance)
:
    stencil(stencil),
    dimensions(dimensions),
    maxTerms(Polynomial::nTerms(dimensions)),
    tolerance(tolerance)
{}

template<class Polynomial>
autoPtr<scalarRectangularMatrix> Foam::AdaptivePolynomial<Polynomial>::matrix() const
{
    List<label> fittableTerms(0, label(0));
    for (int terms=0; terms<maxTerms; terms++)
    {
        scalarRectangularMatrix B(stencil.size(), fittableTerms.size()+1, scalar(0));
        for (int i=0; i<B.n(); i++)
        {
            scalar coefficients[maxTerms]; 
            Polynomial::addCoeffs(coefficients, stencil[i], 1, dimensions);

            int col = 0;
            for (int j=0; j<maxTerms; j++)
            {
                if (containsEntry(fittableTerms, j) || j == terms)
                {
                    B[i][col++] = coefficients[j];
                }
            }
        }
        if (fullRank(B))
        {
            fittableTerms.append(terms);
        }
    }
    
    scalarRectangularMatrix* B = new scalarRectangularMatrix(stencil.size(), fittableTerms.size(), scalar(0));
    for (int i=0; i<B->n(); i++)
    {
        scalar coefficients[maxTerms]; 
        Polynomial::addCoeffs(coefficients, stencil[i], 1, dimensions);
        int col = 0;
        for (int j=0; j<maxTerms; j++)
        {
            if (containsEntry(fittableTerms, j))
            {
                (*B)[i][col++] = coefficients[j];
            }
        }
    }

    return autoPtr<scalarRectangularMatrix>(B);
}

template<class Polynomial>
bool Foam::AdaptivePolynomial<Polynomial>::containsEntry(const List<label> fittableTerms, label term) const
{
    return findIndex(fittableTerms, term) != -1;
}

template<class Polynomial>
bool Foam::AdaptivePolynomial<Polynomial>::fullRank(const scalarRectangularMatrix& B) const
{
    SVD svd(B, SMALL);
    return min(svd.S()) > tolerance;
}

