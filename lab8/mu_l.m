function [mu] = mu_l(z,K)
    if (z >=127 + K)
        mu = 1;
    elseif ((z>127) && (z<127+K))
        mu = (z - 127) / K;
    else 
        mu = 0;
    end
end