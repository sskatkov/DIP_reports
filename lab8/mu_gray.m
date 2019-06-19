function mu = mu_gray( z, K, mid )
    a = mid;
    if (z < a) && (z >= a-K)
        mu = 1 - (a-z)/K;
    elseif (z >= a) && (z <= a+K)
        mu = 1 - (z-a)/K;
    else
        mu = 0;
    end
end

