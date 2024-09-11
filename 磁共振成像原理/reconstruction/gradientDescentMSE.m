function rho = gradientDescent(b,E,maxit)
    siz(1) = size(b,1);
    siz(2) = size(b,2);
    m = E'*b;
    r = E'*b - E'*(E*m);
    rprev = r;
    for it = 1:maxit
        aux = E'*(E*r);
        alpha = (r(:)'*r(:)) / (r(:)'*aux(:));
        m = m + alpha*r;
        r = r - alpha*(E'*(E*r));
        err = sum(abs(r(:)));
        rprev = r;
        rho = m;
        if err < 5
            it
            rho = m;
            break
        end

    end
end