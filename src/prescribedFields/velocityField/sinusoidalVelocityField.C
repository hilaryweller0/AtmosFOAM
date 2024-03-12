#include "sinusoidalVelocityField.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(sinusoidalVelocityField, 0);
addToRunTimeSelectionTable
(
    velocityField,
    sinusoidalVelocityField,
    dict
);

sinusoidalVelocityField::sinusoidalVelocityField
(
    const dictionary& dict
)
:
    a_x(dict.lookupOrDefault<scalar>("a_x", scalar(0))),
    a_y(dict.lookupOrDefault<scalar>("a_y", scalar(0))),
    a_z(dict.lookupOrDefault<scalar>("a_z", scalar(0))),
    
    b_x(dict.lookupOrDefault<scalar>("b_x", scalar(1))),
    b_y(dict.lookupOrDefault<scalar>("b_y", scalar(1))),
    b_z(dict.lookupOrDefault<scalar>("b_z", scalar(1))),
    
    c_x(dict.lookupOrDefault<scalar>("c_x", scalar(0))),
    c_y(dict.lookupOrDefault<scalar>("c_y", scalar(0))),
    c_z(dict.lookupOrDefault<scalar>("c_z", scalar(0))),
    
    d_x(dict.lookupOrDefault<scalar>("d_x", scalar(1))),
    d_y(dict.lookupOrDefault<scalar>("d_y", scalar(1))),
    d_z(dict.lookupOrDefault<scalar>("d_z", scalar(1))),
    
    e_x(dict.lookupOrDefault<scalar>("e_x", scalar(0))),
    e_y(dict.lookupOrDefault<scalar>("e_y", scalar(0))),
    e_z(dict.lookupOrDefault<scalar>("e_z", scalar(0))),
    
    f_x(dict.lookupOrDefault<scalar>("f_x", scalar(1))),
    f_y(dict.lookupOrDefault<scalar>("f_y", scalar(1))),
    f_z(dict.lookupOrDefault<scalar>("f_z", scalar(1))),
    
    k_x(dict.lookupOrDefault<scalar>("k_x", scalar(0))),
    k_y(dict.lookupOrDefault<scalar>("k_y", scalar(0))),
    k_z(dict.lookupOrDefault<scalar>("k_z", scalar(0))),
    
    l_x(dict.lookupOrDefault<scalar>("l_x", scalar(0))),
    l_y(dict.lookupOrDefault<scalar>("l_y", scalar(0))),
    l_z(dict.lookupOrDefault<scalar>("l_z", scalar(0))),
    
    m_x(dict.lookupOrDefault<scalar>("m_x", scalar(0))),
    m_y(dict.lookupOrDefault<scalar>("m_y", scalar(0))),
    m_z(dict.lookupOrDefault<scalar>("m_z", scalar(0))),
    
    L_x(dict.lookupOrDefault<scalar>("L_x", scalar(1))),
    L_y(dict.lookupOrDefault<scalar>("L_y", scalar(1))),
    L_z(dict.lookupOrDefault<scalar>("L_z", scalar(1))),
    
    phi_x(dict.lookupOrDefault<scalar>("phi_x", scalar(0))),
    phi_y(dict.lookupOrDefault<scalar>("phi_y", scalar(0))),
    phi_z(dict.lookupOrDefault<scalar>("phi_z", scalar(0))),
    
    psi_x(dict.lookupOrDefault<scalar>("psi_x", scalar(0))),
    psi_y(dict.lookupOrDefault<scalar>("psi_y", scalar(0))),
    psi_z(dict.lookupOrDefault<scalar>("psi_z", scalar(0))),
    
    chi_x(dict.lookupOrDefault<scalar>("chi_x", scalar(0))),
    chi_y(dict.lookupOrDefault<scalar>("chi_y", scalar(0))),
    chi_z(dict.lookupOrDefault<scalar>("chi_z", scalar(0)))
{};

vector sinusoidalVelocityField::velocityAt
(
    const point& p,
    const Time& t
) const
{
    const dimensionedScalar T = t.endTime();

    scalar u = (a_x + b_x*Foam::cos(2*M_PI*k_x*(p.x()-phi_x)/L_x))
             * (a_y + b_y*Foam::cos(2*M_PI*k_y*(p.y()-phi_y)/L_y))
             * (a_z + b_z*Foam::cos(2*M_PI*k_z*(p.z()-phi_z)/L_z));
    scalar v = (c_x + d_x*Foam::cos(2*M_PI*l_x*(p.x()-psi_x)/L_x))
             * (c_y + d_y*Foam::cos(2*M_PI*l_y*(p.y()-psi_y)/L_y))
             * (c_z + d_z*Foam::cos(2*M_PI*l_z*(p.z()-psi_z)/L_z));
    scalar w = (e_x + f_x*Foam::cos(2*M_PI*m_x*(p.x()-chi_x)/L_x))
             * (e_y + f_y*Foam::cos(2*M_PI*m_y*(p.y()-chi_y)/L_y))
             * (e_z + f_z*Foam::cos(2*M_PI*m_z*(p.z()-chi_z)/L_z));

    return vector(u,v,w);
}

}
