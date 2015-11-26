#include "Mountain.H"
#include "fvCFD.H"

SchaerFineMountain::SchaerFineMountain(const scalar lambda) : lambda(lambda) {};

scalar SchaerFineMountain::heightAt(const point& p) const
{
        return sqr(Foam::cos(M_PI * p.x() / lambda));
}
