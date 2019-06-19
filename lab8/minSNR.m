function [rsnr] = minSNR(D_0, I, J, U, V)
    H = exp(-(U.^2+V.^2)/(2 * D_0^2));
    F = fft2(J, size(H, 1), size(H, 2));
    g = real(ifft2(H .* F));
    g = (g(1:size(J, 1), 1:size(J, 2)));
    rsnr = snr(double(I), double(g) - double(I));
end
