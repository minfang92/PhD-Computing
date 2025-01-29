using Revise
includet("../src/FHANK.jl")

#Step1: Choose some parameters
p = Param()
p.Πz
pw = 1.0
w = 1.2 # necessary to start with low wage...
Λ = p.β
U = 0
q=1.0
typeof(U)

@time Firm_VFI(p, U, pw, w, Λ, q)




p.ξbar

p.vf
p.vfA
p.vfNA
p.ξstar
p.kpolicyA
p.kpolicyNA

#Step2: Setup reset function
function vfi_reset(p::NamedTuple)
        for ik in 1:p.Nk
            for iz in 1:p.Nz
               p.vf[ik,iz] = ik * iz / (ik+iz)
               p.kpolicy[ik,iz] = ik/10
            end
        end
end

#Step3.1: Test laborpolicy
lll = laborpolicy(p, S, pw, w)
profit = profit_func(p, S, pw, w)


#Step3.2: Test capitalpolicy
iz = 3
ik = 2
EV = expectation(p, p.vf, p.Πz)
EV_func = interpolate(p.kgrid, EV[:,iz],SteffenMonotonicInterpolation())
l = p.lpolicy[ik,iz]
zprod = exp(p.zgrid[iz])
k = p.kgrid[ik]
y = zprod * k^p.α * l^p.ν
πl = pw * y - w * l

#Step4: Test Expectation and Bellman
vfi_reset(p)
p.vf
EV = expectation(p, p.vf, p.Πz)
bellman(p, pw, w, Λ, EV, profit)
p.vf

#Step5: Test Howard
p.kpolicy
howard(p, S, pw, w, Λ)
p.vf

#Step6: Test VFI
vfi_reset(p)
p.vf

p.vf
p.kpolicy

plot(p.kgrid,p.kpolicy[:,1])
