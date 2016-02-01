#ifndef WITCH_OF_AGNESI_MOUNTAIN_H
#define WITCH_OF_AGNESI_MOUNTAIN_H

#include "fvCFD.H"
#include "Mountain.H"

class WitchOfAgnesi3DMountain : public Mountain
{
    public:
    TypeName("WitchOfAgnesi3DMountain");
    WitchOfAgnesi3DMountain(const IOdictionary& dict);
    scalar heightAt(const point& p) const;
    scalar gradientAt(const scalar x) const;
    scalar timeToCross(const scalar u0, const scalar H) const;
    scalar halfWidth() const;
    private:
    const scalar a;
    const scalar h0;
};

#endif
