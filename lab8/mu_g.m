function [mu] = mu_g(z,K)
    if ((z >= 127-K) && (z <= 127))
        mu = (z - 127 + K) / K;
    elseif ((z > 127) && (z <= 127 + K))
        mu = -(z - 127 - K) / K;
    else
        mu = 0;
    end
end