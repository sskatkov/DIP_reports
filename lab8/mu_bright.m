function mu = mu_bright( z, K, mid )
    a = mid;
    if (z <= a) && (z >= a-K)
        mu = 1 - (a-z)/K;
    elseif z > a
        mu = 1;
    else
        mu = 0;
    end
end

