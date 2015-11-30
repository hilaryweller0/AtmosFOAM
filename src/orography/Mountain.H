#ifndef MOUNTAIN_H
#define MOUNTAIN_H

#include "fvCFD.H"

class AbstractMountain
{
    public:
    virtual scalar heightAt(const point& p) const = 0;
    virtual ~AbstractMountain() {};
};

class SchaerFineMountain : public AbstractMountain
{
    public:
    SchaerFineMountain(const scalar lambda);
    scalar heightAt(const point& p) const;
    private:
    const scalar lambda;
};

class SchaerSmoothCosMountain : public AbstractMountain
{
    public:
    SchaerSmoothCosMountain(const scalar a, const scalar h0);
    scalar heightAt(const point& p) const;
    private:
    const scalar a;
    const scalar h0;
};

class SchaerSmoothExpMountain : public AbstractMountain
{
    public:
    SchaerSmoothExpMountain(const scalar a, const scalar h0);
    scalar heightAt(const point& p) const;
    private:
    const scalar a;
    const scalar h0;
};

class Mountain : public AbstractMountain
{
    public:
    declareRunTimeSelectionTable
    (
        autoPtr,
        Mountain,
        dict,
        (
            const IOdictionary& dict
        ),
        (dict)
    );

    static autoPtr<Mountain> New(const IOdictionary& dict);
    virtual ~Mountain() {};

    virtual scalar gradientAt(const scalar x) const = 0;
    virtual scalar timeToCross(const scalar u0, const scalar H) const = 0;
    virtual scalar halfWidth() const = 0;
};

#endif
