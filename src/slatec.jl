module SLATEC

export DRC, DRD, DRF, DRJ

const D1MACH1 = floatmin(Float64)
const D1MACH2 = floatmax(Float64)
const D1MACH3 = eps(Float64)/2
const D1MACH4 = eps(Float64)
const D1MACH5 = log10(2.)

#***BEGIN PROLOGUE  DRF
#***PURPOSE  Compute the incomplete or complete elliptic integral of the
#            1st kind.  For X, Y, and Z non-negative and at most one of
#            them zero, RF(X,Y,Z) = Integral from zero to infinity of
#                                -1/2     -1/2     -1/2
#                      (1/2)(t+X)    (t+Y)    (t+Z)    dt.
#            If X, Y or Z is zero, the integral is complete.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RF-S, DRF-D)
#***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
#             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE FIRST KIND,
#             TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

function DRF(X, Y, Z)

    ERRTOL = (4.0*D1MACH3)^(1.0/6.0)
    LOLIM  = 5.0 * D1MACH1
    UPLIM  = D1MACH2/5.0
    C1 = 1.0/24.0
    C2 = 3.0/44.0
    C3 = 1.0/14.0

    ans = 0.0
    if min(X,Y,Z) < 0.0
        return ans, 1
    end

    if max(X,Y,Z) > UPLIM
        return ans, 3
    end

    if min(X+Y,X+Z,Y+Z) < LOLIM
        return ans, 2
    end

    XN = X
    YN = Y
    ZN = Z
    MU = 0.
    XNDEV = 0.
    YNDEV = 0.
    ZNDEV = 0.

    while true
        MU = (XN+YN+ZN)/3.0
        XNDEV = 2.0 - (MU+XN)/MU
        YNDEV = 2.0 - (MU+YN)/MU
        ZNDEV = 2.0 - (MU+ZN)/MU
        EPSLON = max(abs(XNDEV),abs(YNDEV),abs(ZNDEV))
        if (EPSLON < ERRTOL) break end
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
        XN = (XN+LAMDA)*0.250
        YN = (YN+LAMDA)*0.250
        ZN = (ZN+LAMDA)*0.250
    end

    E2 = XNDEV*YNDEV - ZNDEV*ZNDEV
    E3 = XNDEV*YNDEV*ZNDEV
    S  = 1.0 + (C1*E2-0.10-C2*E3)*E2 + C3*E3
    ans = S/sqrt(MU)

    return ans, 0
end

#***BEGIN PROLOGUE  DRD
#***PURPOSE  Compute the incomplete or complete elliptic integral of
#            the 2nd kind. For X and Y nonnegative, X+Y and Z positive,
#            DRD(X,Y,Z) = Integral from zero to infinity of
#                                -1/2     -1/2     -3/2
#                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.
#            If X or Y is zero, the integral is complete.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RD-S, DRD-D)
#***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
#             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE SECOND KIND,
#             TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

function DRD(X, Y, Z)

    ERRTOL = (D1MACH3/3.0)^(1.0/6.0)
    LOLIM  = 2.0/(D1MACH2)^(2.0/3.0)
    TUPLIM = D1MACH1^(1.0E0/3.0E0)
    TUPLIM = (0.10*ERRTOL)^(1.0E0/3.0E0)/TUPLIM
    UPLIM  = TUPLIM^2.0
    C1 = 3.0/14.0
    C2 = 1.0/6.0
    C3 = 9.0/22.0
    C4 = 3.0/26.0

    ans = 0.0
    if min(X,Y) < 0.0
        return ans, 1
    end

    if max(X,Y,Z) > UPLIM
        return ans, 3
    end

    if min(X+Y,Z) < LOLIM
        return ans, 2
    end

    XN = X
    YN = Y
    ZN = Z
    SIGMA = 0.0
    POWER4 = 1.0
    MU = 0.
    XNDEV = 0.
    YNDEV = 0.
    ZNDEV = 0.

    while true
        MU = (XN+YN+3.0*ZN)*0.20
        XNDEV = (MU-XN)/MU
        YNDEV = (MU-YN)/MU
        ZNDEV = (MU-ZN)/MU
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))
        if (EPSLON < ERRTOL) break end
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
        SIGMA = SIGMA + POWER4/(ZNROOT*(ZN+LAMDA))
        POWER4 = POWER4*0.250
        XN = (XN+LAMDA)*0.250
        YN = (YN+LAMDA)*0.250
        ZN = (ZN+LAMDA)*0.250
    end

    EA = XNDEV*YNDEV
    EB = ZNDEV*ZNDEV
    EC = EA - EB
    ED = EA - 6.0*EB
    EF = ED + EC + EC
    S1 = ED*(-C1+0.250*C3*ED-1.50*C4*ZNDEV*EF)
    S2 = ZNDEV*(C2*EF+ZNDEV*(-C3*EC+ZNDEV*C4*EA))
    ans = 3.0*SIGMA + POWER4*(1.0+S1+S2)/(MU*sqrt(MU))

    return ans, 0
end

#***BEGIN PROLOGUE  DRC
#***PURPOSE  Calculate a double precision approximation to
#             DRC(X,Y) = Integral from zero to infinity of
#                              -1/2     -1
#                    (1/2)(t+X)    (t+Y)  dt,
#            where X is nonnegative and Y is positive.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RC-S, DRC-D)
#***KEYWORDS  DUPLICATION THEOREM, ELEMENTARY FUNCTIONS,
#             ELLIPTIC INTEGRAL, TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

function DRC(X, Y)

    ERRTOL = (D1MACH3/16.0)^(1.0/6.0)
    LOLIM  = 5.0 * D1MACH1
    UPLIM  = D1MACH2 / 5.0
    C1 = 1.0/7.0
    C2 = 9.0/22.0

    ans = 0.0
    if X < 0.0 || Y <= 0.0
        return ans, 1
    end

    if max(X,Y) > UPLIM
        return ans, 3
    end

    if X+Y < LOLIM
        return ans, 2
    end

    XN = X
    YN = Y
    MU = 0.
    SN = 0.

    while true
        MU = (XN+YN+YN)/3.0
        SN = (YN+MU)/MU - 2.0
        if abs(SN) < ERRTOL break end
        LAMDA = 2.0*sqrt(XN)*sqrt(YN) + YN
        XN = (XN+LAMDA)*0.250
        YN = (YN+LAMDA)*0.250
    end

    S = SN*SN*(0.30+SN*(C1+SN*(0.3750+SN*C2)))
    ans = (1.0+S)/sqrt(MU)

    return ans,0
end

#***BEGIN PROLOGUE  DRJ
#***PURPOSE  Compute the incomplete or complete (X or Y or Z is zero)
#            elliptic integral of the 3rd kind.  For X, Y, and Z non-
#            negative, at most one of them zero, and P positive,
#             RJ(X,Y,Z,P) = Integral from zero to infinity of
#                              -1/2     -1/2     -1/2     -1
#                    (3/2)(t+X)    (t+Y)    (t+Z)    (t+P)  dt.
#***LIBRARY   SLATEC
#***CATEGORY  C14
#***TYPE      DOUBLE PRECISION (RJ-S, DRJ-D)
#***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
#             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE THIRD KIND,
#             TAYLOR SERIES
#***AUTHOR  Carlson, B. C.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Notis, E. M.
#             Ames Laboratory-DOE
#             Iowa State University
#             Ames, IA  50011
#           Pexton, R. L.
#             Lawrence Livermore National Laboratory
#             Livermore, CA  94550

function DRJ(X, Y, Z, P)

    ERRTOL = (D1MACH3/3.0)^(1.0/6.0)
    LOLIM  = (5.0 * D1MACH1)^(1.0/3.0)
    UPLIM  = 0.30*( D1MACH2 / 5.0)^(1.0/3.0)

    C1 = 3.0/14.0
    C2 = 1.0/3.0
    C3 = 3.0/22.0
    C4 = 3.0/26.0

    ans = 0.0
    if min(X,Y,Z) < 0.0
        return ans, 1
    end

    if max(X,Y,Z,P) > UPLIM
        return ans, 3
    end

    if min(X+Y,X+Z,Y+Z,P) < LOLIM
        return ans, 2
    end

    IER = 0
    XN = X
    YN = Y
    ZN = Z
    PN = P
    SIGMA = 0.0
    POWER4 = 1.0
    MU = 0.
    XNDEV = 0.
    YNDEV = 0.
    ZNDEV = 0.
    PNDEV = 0.

    while true
        MU = (XN+YN+ZN+PN+PN)*0.20
        XNDEV = (MU-XN)/MU
        YNDEV = (MU-YN)/MU
        ZNDEV = (MU-ZN)/MU
        PNDEV = (MU-PN)/MU
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV), abs(PNDEV))
        if EPSLON < ERRTOL break end
        XNROOT =  sqrt(XN)
        YNROOT =  sqrt(YN)
        ZNROOT =  sqrt(ZN)
        LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
        ALFA = PN*(XNROOT+YNROOT+ZNROOT) + XNROOT*YNROOT*ZNROOT
        ALFA = ALFA*ALFA
        BETA = PN*(PN+LAMDA)*(PN+LAMDA)
        drc,IER = DRC(ALFA,BETA)
        SIGMA = SIGMA + POWER4*drc
        POWER4 = POWER4*0.250
        XN = (XN+LAMDA)*0.250
        YN = (YN+LAMDA)*0.250
        ZN = (ZN+LAMDA)*0.250
        PN = (PN+LAMDA)*0.250
    end

    EA = XNDEV*(YNDEV+ZNDEV) + YNDEV*ZNDEV
    EB = XNDEV*YNDEV*ZNDEV
    EC = PNDEV*PNDEV
    E2 = EA - 3.0*EC
    E3 = EB + 2.0*PNDEV*(EA-EC)
    S1 = 1.0 + E2*(-C1+0.750*C3*E2-1.50*C4*E3)
    S2 = EB*(0.50*C2+PNDEV*(-C3-C3+PNDEV*C4))
    S3 = PNDEV*EA*(C2-PNDEV*C3) - C2*PNDEV*EC
    ans = 3.0*SIGMA + POWER4*(S1+S2+S3)/(MU* sqrt(MU))

    return ans, IER
end

end # module
