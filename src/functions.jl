function B_mixed(S; mu = -0.1, sigma = 0.05)
    B = rand(Normal(mu, sigma), S, S)
    B[diagind(B)] .= -1
    B
end

function glv(u, p, t)
    r, A, K = p
    r .* (1 .+ A * u ./ K) .* u
end
