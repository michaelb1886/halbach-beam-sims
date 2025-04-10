using LinearAlgebra

function field_from_mag_dipole(r,m)
    """
    returns the field from a magnetic dipole moment m
    at the position r 
    SI units (m: [J/T]=[Am^2], r: [m])
    """

    rn = norm(r)
    1.0e-7*(3*r.*dot(m,r)./(rn^5) .- m./(rn^3))
end

