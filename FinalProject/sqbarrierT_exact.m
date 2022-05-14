function T = sqbarrierT_exact(kb,eta)
    xib = sqrt(eta-kb.^2);
    S = sinh(2*xib);
    T = 1./(1 + S.^2.*(kb./xib + xib./kb).^2/4);
