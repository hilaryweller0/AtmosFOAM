#ifndef MOUNTAIN_H
#define MOUNTAIN_H

#include "fvCFD.H"

class Mountain
{
    public:
    virtual scalar heightAt(const scalar x) const = 0;
    virtual ~Mountain() {};
};

class SchaerCosFineMountain : public Mountain
{
    public:
    SchaerCosFineMountain(const scalar lambda);
    scalar heightAt(const scalar x) const;
    private:
    const scalar lambda;
};

class SchaerSmoothMountain : public Mountain
{
    public:
    SchaerSmoothMountain(const scalar a, const scalar h0);
    scalar heightAt(const scalar x) const;
    private:
    const scalar a;
    const scalar h0;
};

class SchaerCosMountain : public Mountain
{
    public:
    SchaerCosMountain(const IOdictionary& dict);
    SchaerCosMountain(const scalar a, const scalar h0, const scalar lambda);
    scalar heightAt(const scalar x) const;
    scalar gradientAt(const scalar x) const;
    ~SchaerCosMountain();
    private:
    const scalar a;
    const scalar h0;
    const scalar lambda;
    SchaerSmoothMountain *smooth;
    SchaerCosFineMountain *fine;
};

#endif
