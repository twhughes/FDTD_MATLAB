function gs = pulse(trange, fmax, deltat)

    tau = 0.5/fmax;
    t0 = tau*6;   %#recommended starting time for source
    gs = [];
    for t = trange
        gs = [gs, g(t,t0,tau, deltat)];
    end
end
    
function src = g(t, t0, tau, deltat)
    src = exp(-1*((t-t0+deltat)^2)/(tau^2));
end
