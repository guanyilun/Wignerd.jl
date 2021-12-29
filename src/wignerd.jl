module wignerd

using Tullio

alpha(l, s1, s2) = (l ≤ abs(s1) || l ≤ abs(s2)) ? 0. : sqrt((l^2-s1^2)*(l^2-s2^2))/l

function wigd_init!(s1, s2, cosθ, out)
    s12sign = (mod(s1+s2,2) == 1) ? -1 : 1
    A = 1.  # prefactor

    if abs(s1) > abs(s2)
        A *= s12sign
        s1, s2 = s2, s1
    end

    if s2 < 0
        A *= s12sign
        s1, s2 = -s1, -s2
    end

    abs_s1 = abs(s1);

    for i = 1:s2-abs_s1
        A *= sqrt((s2+abs_s1+i)/i)
    end

    @. out = A*((1+cosθ)/2)^((s2+s1)/2) * ((1-cosθ)/2)^((s2-s1)/2)

    s2
end

function wigd_rec!(l, s1, s2, cosθ, wigd_hi, wigd_lo)
    alpha_hi = alpha(l+1, s1, s2)
    alpha_lo = alpha(l, s1, s2)
    beta = (s1==0 || s2==0) ? 0. : (s1*s2/(l*(l+1)))
    @inbounds @simd for i in eachindex(cosθ)
        x = (2*l+1)*(cosθ[i]-beta) * wigd_hi[i] - alpha_lo*wigd_lo[i]
        wigd_lo[i] = wigd_hi[i]
        wigd_hi[i] = x / alpha_hi
    end
end

function cf_from_cl(s1, s2, lmax, cl::AbstractArray{T,1}, cosθ) where T
    wigd_lo = zero(cosθ)
    wigd_hi = zero(cosθ)
    cf      = zero(cosθ)

    l = wigd_init!(s1, s2, cosθ, wigd_hi)
    if l ≤ lmax; cf .= cl[l+1] .* wigd_hi end

    while l < lmax
        wigd_rec!(l, s1, s2, cosθ, wigd_hi, wigd_lo)
        l += 1
        cf .+= cl[l+1] .* wigd_hi
    end
    cf
end

# processing multiple ps at the same time, for nspec=2 it's faster by ~10%,
# for nspec=5 it's faster by 50%. The improvement is marginal because of poor
# memory access pattern when nspec is small.
# input: cl has shape (nell, nspec)
# output: cf has shape (ntheta, nspec)
function cf_from_cl(s1, s2, lmax, cl::AbstractArray{T,2}, cosθ) where T
    wigd_lo = zero(cosθ)
    wigd_hi = zero(cosθ)
    cf      = zeros(T, length(cosθ), size(cl,2))

    l₀ = wigd_init!(s1, s2, cosθ, wigd_hi)
    if l₀ ≤ lmax; (v=view(cl,l₀+1,:); @tullio cf[j,i] = v[i] * wigd_hi[j]) end

    for l = l₀:lmax-1
        wigd_rec!(l, s1, s2, cosθ, wigd_hi, wigd_lo)
        v=view(cl,l+2,:); @tullio cf[j,i] += v[i] * wigd_hi[j]
    end
    cf
end

function cl_from_cf(s1, s2, lmax, cf::AbstractArray{T,1}, cosθ, weights) where T
    wigd_lo = zero(cosθ)
    wigd_hi = zero(cosθ)
    cl = zeros(Float64, lmax+1)

    l = wigd_init!(s1, s2, cosθ, wigd_hi)
    wigd_hi .*=  weights

    if l ≤ lmax; cl[l+1] = (@tullio c=cf[i]*wigd_hi[i]) end

    while l < lmax
        wigd_rec!(l, s1, s2, cosθ, wigd_hi, wigd_lo)
        l += 1
        cl[l+1] = (@tullio c=cf[i]*wigd_hi[i])
    end
    cl
end

# processing multiple cf at the same time, for nspec=2 it's faster by 30%,
# for nspec=5 it's ~ 3 times faster, compared to repeated calls.
# input: cf has shape (ntheta, nspec)
# output: cl has shape (nell, nspec)
function cl_from_cf(s1, s2, lmax, cf::AbstractArray{T,2}, cosθ, weights) where T
    wigd_lo = zero(cosθ)
    wigd_hi = zero(cosθ)
    cl = zeros(Float64, lmax+1, size(cf,2))

    l = wigd_init!(s1, s2, cosθ, wigd_hi)
    wigd_hi .*=  weights

    if l ≤ lmax; (v = view(cl,l+1,:); @tullio v[j] = cf[i,j]*wigd_hi[i]) end

    while l < lmax
        wigd_rec!(l, s1, s2, cosθ, wigd_hi, wigd_lo)
        l += 1
        v = view(cl,l+1,:); @tullio v[j] = cf[i,j]*wigd_hi[i]
    end
    cl
end

end # module
