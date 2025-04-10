using LinearAlgebra, Rotations, Plots

function field_from_mag_dipole(r,m,r0=[0.0,0.0,0.0])
    """
    returns the field from a magnetic dipole moment m
    at the position r 
    SI units (m: [J/T]=[Am^2], r: [m])
    """
    r = r .- r0
    rn = norm(r)
    1.0e-7*(3*r.*dot(m,r)./(rn^5) .- m./(rn^3))
end

function field_from_n_mag_dipole(r, m, r0)
    """
    returns the field from a magnetic dipole moment m
    at the position r 
    SI units (m: [J/T]=[Am^2], r: [m])
    """
    b = zeros(3)
    for (i,mi) in enumerate(m)
        b .+= field_from_mag_dipole(r,mi,r0[i])
    end
    return b
end

spacing = 250e-2
n = 5
angle_shift = pi/2
Br = 1.3
rm,hm  = 100e-3,50e-3
V = pi*rm^2*hm
mu0 = 4*pi*1.0e-7
m0 = [0,1.0e-7*Br*V/mu0,0]
angle0 = 0.0

positions = [LinRange(-(n-1)*spacing/2,(n-1)*spacing/2,n) zeros(n) zeros(n)]
moments = [RotZ(angle_shift*n+angle0)*m0 for i in 0:n-1]

# Ngrid=128
# b = zeros(Ngrid,Ngrid,3)
# x = LinRange(-15,15,Ngrid)
# y = LinRange(-15,15,Ngrid)

# for (xi,x) in enumerate(x)
#     for (yi,y) in enumerate(y)
#         b[xi,yi,:] = field_from_n_mag_dipole([x,y,0], moments, positions)
#     end
# end

# # xgrid = cat([x for i in 1:Ngrid]...;dims=2)
# # ygrid = cat([y for i in 1:Ngrid]...;dims=2)'
# # 
# bmag = dropdims(mapslices(norm, b, dims=3);dims=3)
# p1 = heatmap(x,y,log10.(bmag),
#         xlabel="x [m]", ylabel="y [m]",
#         aspect_ratio=1,
#         xlims=(-15,15), ylims=(-15,15),
#         size=(800,800),
#         colorbar=true,
#         grid=false,
#         framestyle=:box,)

nθ = 128
θ = LinRange(0,2π,nθ)
p = plot()
bp = zeros(nθ,3)
bpm = zeros(nθ)
for r in [1,5,10,50] #LinRange(1,10,5)
# r = 1
    for θi in eachindex(θ)
        x = r*cos(θ[θi])
        y = r*sin(θ[θi])
        bp[θi,:] = field_from_n_mag_dipole([x,y,0], moments, positions)
    end
    bpm[:] = log10.(dropdims(mapslices(norm, bp;dims=2);dims=2))
    plot!(p, θ, (bpm).-minimum(bpm), proj = :polar)
end
display(p)


# new figure