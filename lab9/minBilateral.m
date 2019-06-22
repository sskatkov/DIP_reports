function [rsnr]  = minBilateral(I, J,sigma_d, sigma_r) 
    if (sigma_d > 0 && sigma_r)
        Out = bilateral_filter1(J,fix(3*sigma_d+1),sigma_d,sigma_r);
        rsnr = snr(double(I), double(Out) - double(I));
    else
        rsnr = 0;
    end
end