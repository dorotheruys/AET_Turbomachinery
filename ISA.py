import math as m
g0 = float(9.80665)
R = float(287)
gam = float(1.4)
Tb1 = float(288.15)
pb1 = float(101325)

def calc(Tb,a,hl,pb):
    T = Tb + a * hl
    p = pb * ((T / Tb) ** ((-1 * g0) / (a * R)))
    rho = (p) / (R * T)
    sof = m.sqrt(R * gam * T)
    return(T, p, rho, sof)

def ISA(h):
    if h <= 11000:
        hl = h
        a = float(-0.0065)
        Tb = float(Tb1)
        pb = float(pb1)
        print("You've reached the Troposhere.")
        res = calc(Tb,a,hl,pb)
        return(res)