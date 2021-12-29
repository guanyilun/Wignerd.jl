module wignerd

function dot(a, b)
    s = 0
    @inbounds @simd for i in eachindex(a)
        s += a[i] * b[i]
    end
    s
end
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

function wigd_cf_from_cl(s1, s2, lmax, cl, cosθ)
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
end

function wigd_cl_from_cf(s1, s2, lmax, cf, cosθ, weights)
    wigd_lo = zero(cosθ)
    wigd_hi = zero(cosθ)
    cl = zeros(Float64, lmax+1)
    
    l = wigd_init!(s1, s2, cosθ, wigd_hi)
    wigd_hi .*=  weights
    
    if l ≤ lmax; cl[l+1] = dot(cf, wigd_hi) end

    while l < lmax
        wigd_rec!(l, s1, s2, cosθ, wigd_hi, wigd_lo)
        l += 1
        cl[l+1] = dot(cf, wigd_hi)
    end
    cl
end

end # module
