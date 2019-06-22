function R = fuzzy(z, M)
    L = 256;
    d = z' - z(5); 
    nu_ = @(d) (d <= M) * acot(M+L-1)*(M-d);
    nu = @(d) (d >= -M) * acot(M+L-1)*(M+d);
    l_ = min([nu_(d(2)) nu_(d(4)) nu_(d(6)) nu_(d(8))]);
    l = min([nu(d(2)) nu(d(4)) nu(d(6)) nu(d(8))]);
    if l_ == 0 && l == 0 
        R = z(5);
        return; 
    end
    nu_ = @(v) (v <= 0) * -acot(L)*v;
    nu = @(v) (v >= 0) * acot(L)*v;
    Q = @(v) max([min([l_ nu_(v)]) min([l nu(v)])]);
    num = 0;
    den = 0;
    for v = -L:L 
        num = num + v*Q(v);
        den = den + Q(v);
    end
    v0 = num / den; 
    R = z(5) + v0;
end
