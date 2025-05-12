function eqm = Model_MarketShare(p, np, func, eqm)

    eqm.s(1) = (1 - p.gamma) * eqm.s(2) * eqm.mu_hat(1,2) + (1 - p.gamma) * eqm.s(1) - (1 - p.gamma) * eqm.s(1) * eqm.mu_hat(2,1) + p.gamma * p.m(1);
    eqm.s(2) = (1 - p.gamma) * eqm.s(1) * eqm.mu_hat(2,1) + (1 - p.gamma) * eqm.s(2) - (1 - p.gamma) * eqm.s(2) * eqm.mu_hat(1,2) + p.gamma * p.m(2);


    flag = eqm.s(1) < 0 || eqm.s(1) > 1 || eqm.s(2) < 0 || eqm.s(2) > 1;
    if flag
        disp('Market share cannot be negative or greater than 1.')

        eqm.s(1) = max(0.01, eqm.s(1));
        eqm.s(1) = min(0.99, eqm.s(1));
        eqm.s(2) = 1 - eqm.s(1);
    end



end