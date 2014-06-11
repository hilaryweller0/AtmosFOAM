/* Functions for the different mountain types */

using namespace Foam::constant::mathematical;

Foam::scalar AgnessiWitch(const scalar x, const scalar a, const scalar hm)
{
    return hm/(1 + sqr(x/a));
}

Foam::scalar BottaKlein(const scalar x, const scalar a, const scalar hm)
{
    return hm/Foam::pow(1 + sqr(x/a), 1.5);
}

Foam::scalar ScharExp(const scalar x, const scalar a, const scalar hm)
{
    return hm*Foam::exp(-sqr(x/a));
}

Foam::scalar ScharCos(const scalar x, const scalar a)
{
    return sqr(Foam::cos(M_PI*x/a));
}

Foam::scalar ScharCosSmooth(const scalar x, const scalar a, const scalar hm)
{
    scalar h = 0;
    if (mag(x) < a)
    {
        h = hm*sqr(Foam::cos(0.5*M_PI*x/a));
    }
    return h;
}

Foam::scalar flatMountain(const scalar x, const scalar a)
{
    return scalar(1);
}

