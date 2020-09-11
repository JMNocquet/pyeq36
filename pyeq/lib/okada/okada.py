"""
Wrapper for dc3d0 and dc3d original Okada's routines
This corresponds to the okada.f from edcmp
"""


def okada(llambda, mu, dislocations, xs, ys, zs, lengths, widths, strikes, dips, rakes, xrec, yrec, zrec0):
    """
    c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    c       ns = the really used number of rectangular sources
    c       NSMAX = the upper limit of ns
    c       nrec = the really used number of observation positions
    c       NRECMAX = the upper limit of nrec
    c
    c       lambda, mu = the two Lame constants in Pascal (SI unit)
    c
    c       (xs,ys,zs) = coordinates of the start point of strike
    c       with x = north, y = east, z = downward.
    c       all angles in degree.
    c       (xrec,yrec,zrec0) = cartesian coordinates of observations
    c             Note zrec0 is a fixed constant
    c       disp = 3 displacement components: ux,uy,uz
    c       strain = 6 strain components: exx,eyy,ezz,exy,eyz,ezx
    c       tilt = 2 vertical tilt components: dux/dz, duy/dz
    c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    
    import numpy as np

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#    CHECK FOR VARIABLE TYPES, REQUIRED BY PYTHON
#    ============================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#    LOCAL CONSTANTS
#    ===============
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    degtorad = 1.745329252E-02
    eps = 1.E-6
    
    disp0 = np.zeros(3)
    tilt0 = np.zeros(2)
    strain0 = np.zeros(6)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#    LOCAL WORK SPACES
#    =================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    double precision disp0(3),tilt0(2),strain0(6)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#    PROCESSING
#    ==========
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#    receiver and source independent variables
#
    ALPHA = (llambda + mu) / (llambda + 2.*mu)
    POT3 = 0.0
    POT4 = 0.0
    DISL3 = 0.0
    AL1 = 0.0
    AW2 = 0.0
    Z = -zrec0

    # number of observations
    nrec = xrec.shape[0]
    
    # number of sources
    ns = xs.shape[0]

    #      initialization

    strain = np.zeros((nrec , 6))
    disp = np.zeros((nrec , 3))
    tilt = np.zeros((nrec , 2))

    # loop on observations points
    
    for irec in np.arange(nrec):

        # loop on sources

        for iis in np.arange(ns): 

            st = strikes[iis] * degtorad
            csst = np.cos(st)
            ssst = np.sin(st)
            cs2st = np.cos(2.*st)
            ss2st = np.sin(2.*st)

            di = dips[iis] * degtorad
            csdi = np.cos(di)
            ssdi = np.sin(di)

            ra = rakes[iis] * degtorad
            csra = np.cos(ra)
            ssra = np.sin(ra)

#        transform from Aki's to Okada's system

            X = (xrec[irec] - xs[iis]) * csst + (yrec[irec] - ys[iis]) * ssst
            Y = (xrec[irec] - xs[iis]) * ssst - (yrec[irec] - ys[iis]) * csst
            DEPTH = zs[iis]
            DIP = dips[iis]

############ point source #####################################################
            if (lengths[iis] == 0) and (widths[iis] == 0):

                from pyeq.lib.okada.dc3d0 import dc3d0

                POT1 = dislocations[iis] * csra
                POT2 = dislocations[iis] * ssra

                IRET = 1

                [UX, UY, UZ, UXX, UYX, UZX, UXY, UYY, UZY, UXZ, UYZ, UZZ] = dc3d0(ALPHA, X, Y, Z, DEPTH, DIP, POT1, POT2, POT3, POT4)
#          if(IRET.eq.1)then
#            stop ' There is a problem in Okada subroutine!'
#          endif
            
            else:

########## finite source #######################################################

                from pyeq.lib.okada.dc3d import dc3d

                AL2= lengths[iis]
                AW1=- widths[iis]
                DISL1= dislocations[iis]*csra
                DISL2= dislocations[iis]*ssra 
                if ( lengths[iis] == 0):
                    AL2= widths[iis]*eps
                    DISL1=DISL1/AL2
                    DISL2=DISL2/AL2
                elif ( widths[iis] == 0. ):
                    AW1= lengths[iis]*eps
                    DISL1=DISL1/(-AW1)
                    DISL2=DISL2/(-AW1)

                IRET=1
                [UX, UY, UZ, UXX, UYX, UZX, UXY, UYY, UZY, UXZ, UYZ, UZZ] = \
                        dc3d(ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3)
#          
#if(IRET.eq.1)then
#            stop ' There is a problem in Okada subroutine!'
#          endif
#       endif

#
#        transform from Okada's to Aki's system
#

            disp0[0] = UX * csst + UY * ssst
            disp0[1] = UX * ssst - UY * csst
            disp0[2] = -UZ

            tilt0[0] = -(UXZ * csst + UYZ * ssst)
            tilt0[1] = -(UXZ * ssst - UYZ * csst)

            strain0[0] = UXX * csst * csst + UYY * ssst * ssst + 0.5 * (UXY + UYX) * ss2st
            strain0[1] = UXX * ssst * ssst + UYY * csst * csst - 0.5 * (UXY + UYX) * ss2st
            strain0[2] = UZZ
            strain0[3] =  0.5 * (UXX - UYY) * ss2st - 0.5 * (UXY + UYX) * cs2st
            strain0[4] = -0.5 * (UZX + UXZ) * ssst  + 0.5 * (UYZ + UZY) * csst
            strain0[5] = -0.5 * (UZX + UXZ) * csst  - 0.5 * (UYZ + UZY) * ssst

            disp[ irec, : ] = disp[ irec, : ] + disp0
            tilt[ irec, : ] = tilt[ irec, : ] + tilt0
            strain[ irec, : ] = strain[ irec, : ] + strain0

    return( disp , tilt, strain )
