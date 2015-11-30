#include "Mountain.H"
#include "fvCFD.H"

SchaerFineMountain::SchaerFineMountain(const scalar lambda) : lambda(lambda) {};

scalar SchaerFineMountain::heightAt(const scalar x) const
{
        return sqr(Foam::cos(M_PI * x / lambda));
}
