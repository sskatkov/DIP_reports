function [mu] = mu_d(z,K)
    if (z <= 127-K)
        mu = 1;
    elseif ((z > 127-K) && (z <= 127))
        mu = 1 - (z - 127 + K)/K;
       else mu = 0;
    end
end