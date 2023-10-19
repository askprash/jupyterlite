import numpy as np

def TR(M, gam = 1.4):
    """Returns the total to static temperature ratio given the Mach number """
    return 1 + 0.5*(gam - 1)*M**2

def PR(M, gam = 1.4):
    """Returns the total to static pressure ratio given the Mach number """
    return (1 + 0.5*(gam - 1)*M**2)**(gam/(gam-1))

def D(M, gam = 1.4, R = 287.05):
    """Returns mdot*sqrt(Tt)/(A*Pt); Note R, gamma on RHS"""
    gm1 = gam - 1
    gp1 = gam + 1
    CorrMass_Area = np.sqrt(gam/R)*M*(1+0.5*gm1*M**2)**(-0.5*gp1/gm1)
    return CorrMass_Area

def PRtoM(PR, gam = 1.4):
    gm1 = gam - 1
    M = np.sqrt(2.0/gm1*(PR**(gm1/gam) - 1))
    return M

TSTD = 298.15 #K
PSTD = 101325 #Pa

def Dstd(M, gam = 1.4, R = 287.05):
    """ Returns mdot*sqrt(Tt/Tref)/(A*Pt/Pref); Note R, gamma on RHS
        Returns in units of kg/s """
    gm1 = gam - 1
    gp1 = gam + 1
    CorrMass_Area = np.sqrt(gam/R)*M*(1+0.5*gm1*M**2)**(-0.5*gp1/gm1) * PSTD / np.sqrt(TSTD) #[kg/s /m^2]
    return CorrMass_Area

def Dinvsimple(WcqA):
    M = (-480+np.sqrt(480**2 - 4*(-236.8)*(-4.482 - WcqA)))/(2*-236.8)
    return M

def mdot(Pt, Tt, A, Pamb, gam = 1.4, R = 287.05):
    PR = 0.0
    if Pamb>0.0:
        PR = Pt/Pamb
        M = PRtoM(PR, gam=gam)
    else:
        M = 1.0
    if (M>1.0):
        M = 1.0
    DM = D(M, gam = gam, R = R)
    mdot = DM*A*Pt/np.sqrt(Tt)
    return mdot
