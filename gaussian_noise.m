function res = gaussian_noise(m,P)
    res = chol(P) * randn(size(m,1),1);
end