function [I_contrast] = contrast(I_original, K)
    [M, N] = size(I_original);
    I_contrast = zeros(M,N);
    for i = 1:M
        for j = 1:N
            I_contrast(i,j) = 127 * mu_g(I_original(i,j),K) + 255 * mu_l(I_original(i,j),K);
        end
    end
end