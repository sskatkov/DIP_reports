function [W] =  winer(original, noisy)
    [M, N] = size(original);
    S_eta = var(noisy(:) - original(:));
    S_g = abs(fft2(noisy)/sqrt(M*N)).^2;
    S_f = S_g - S_eta;
    W = S_f ./ (S_f + S_eta);
    W(W < 0) = 0;
end