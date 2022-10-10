"""
2D Elastic Dislocation Library
Uses analytical solution from Segall, P. (2010). Earthquake and volcano deformation. Princeton University Press.
:Note 1: Only dip-slip fault solution implemented so far (08/06/2017)
:Note 2: 0-90 dip means dip to positive axis (right) while -90-0 means dip to the left.
:Note 3: Be careful with signs. For dip slip dislocations, positive slip indicates thrust motion.\
         For a vertical strike-slip fault, positive slip indicates left-lateral motion.


"""

###################################################################
def vss_eq_2d(x,d,s):
###################################################################

    """
    Calculates along strike horizontal surface displacements for a 2D infinitely long vertical strike-slip fault extending from depth d to the surface
    
    :param x     : 1D numpy array of location of points (distance from the fault surface trace)
    :param d     : bottom depth of the fault
    :param s     : slip on the fault (s positive is for left-lateral motion)
    :return ux   : vector of along strike displacements at x, ux = s / np.pi * np.arctan(d/x)
    """
    
    import numpy as np

    ux= s / np.pi * np.arctan(d/x)

    return ux

###################################################################
def vss_buried_eq_2d(x,t,b,s):
###################################################################

    """
    Calculates along strike horizontal surface displacements for a 2D infinitely long vertical strike-slip fault extending from depth t (top) to b (bottom)
    
    :param x     : 1D numpy array of location of points (distance from the fault surface trace)
    :param d     : bottom depth of the fault
    :param s     : slip on the fault (s positive is for left-lateral motion)
    :return ux   : vector of along strike displacements at x, ux =  
    """

    if b <=t:
        raise ValueError('!!! t must be strictly greater than b')


    # computes solution from bottom to surface
    
    ux1= vss_eq_2d(x,b,s)

    # computes solution from top to surface

    ux2= vss_eq_2d(x,t,s)

    # superposition
    
    ux = ux1 - ux2

    return ux

###################################################################
def vss_interseismic_eq_2d(x,d,s):
###################################################################

    """
    Calculates the interseismic displacement along strike horizontal surface displacements for a 2D infinitely long vertical strike-slip fault extending from depth d to the surface
    
    :param x     : 1D numpy array of location of points (distance from the fault surface trace)
    :param d     : bottom depth of the fault
    :param s     : slip on the fault (s positive is for left-lateral motion)
    :return ux   : vector of along strike displacements at x, ux = s / np.pi * np.arctan(x/d)
    """
    
    import numpy as np

    ux= s / np.pi * np.arctan(x/d)

    return ux


###################################################################
def ds_eq_2d(x,d,dip,s):
###################################################################
    """
    Calculates strike-perpendicular horizontal and vertical displacements for a 2D dip slip fault reaching the surface at x=0.
    
    :param x     : 1D numpy array of location of points (distance from the fault surface trace)
    :param d     : bottom depth of the fault
    :param dip   : dip of the fault
    :param s     : slip on the fault (s positive is a reverse fault)
    :return ux,uy: two 1D numpy arrays of horizontal and vertical prediction displacements
    :note: s>0 is for thrust. This is the opposite of Segall's convention
    :ref: Segall, P. (2010). Earthquake and volcano deformation. Princeton University Press. 
    """

    import numpy as np

    
    # switch from Segall's convention
    # in this code, we want s>0 for reverse fault
    
    if dip>0 and dip<90:
        s=-s

    if dip<0 and dip>-90:
        dip=180.+dip
        s=s

    if dip>90:
        s=s


    rdip=np.radians(dip)
    
    cd=np.cos(rdip)
    sd=np.sin(rdip)
    
    pi=np.pi
    
    xd=d/np.tan(rdip)
    
#    xp=2*d/np.sin(2*rdip)
    
    zeta= (x-xd) / d
    
    # horizontal from Segall 2010
    uxs= -s/pi *  ( cd * ( np.arctan( zeta ) - pi/2.*np.sign(x) ) + (sd - zeta * cd ) / (1+zeta**2)  )

    # vertical from Segall 2010
    uys= s/pi *  ( sd * ( np.arctan( zeta ) - pi/2.*np.sign(x) ) + (cd + zeta * sd ) / (1+zeta**2)  )


    # Segall 2010 and Cohen 1999 are equivalent for horizontal displacement
    # However, it looks like Cohen's formula is wrong for vertical displacement
    # horizontal - Cohen, 1999, eq. 22    
    #ux=(s*np.cos(rdip)/np.pi)* \
    #   (\
    #    np.arctan((x-xd)/d) - (x-xp)*d / ((x-xd)**2+d**2) - np.sign(x)*np.pi/2
    #    )
    
    # vertical - Cohen, 1999, eq. 23
    #
    #uy=(s*np.sin(rdip)/np.pi)* \
    #   (\
    #    np.arctan(d/(x-xd)) - x*d / ((x-xd)**2+d**2) - (0.-np.sign(x))*np.pi/2
    #    )
    
    
    return(uxs,uys)


###################################################################
def ds_buried_eq_2d(x,t,b,dip,s,ref_x='surface'):
###################################################################
    """
    Calculates horizontal and vertical displacements for a 2D buried dip slip fault 

    :param x     : numpy array of location of points
    :param t     : top depth of the fault
    :param b     : bottom depth of the fault
    :param dip   : dip of the fault in decimal degrees
    :param s     : slip on the fault (s positive is a reverse fault)
    :param ref_x : reference for x. Can intersection with surface 'surface' (default) or x fault top fault edge 'true_loc'
    :return ux,uy: two 1D numpy arrays of horizontal and vertical prediction displacements
    """

    if b <= t :
        raise ValueError(('!!! ERROR: t must be strictly greater than b: %lf > %lf' % (t,b)))
    
    import numpy as np

    if ref_x == 'true_loc' :
    
        # since the origin for x is the fault top edge
        # we first correct for the fault top
    
        x=x+t/np.tan(np.radians(dip))

    if isinstance(ref_x,float):
        x=x-ref_x+t/np.tan(np.radians(dip))
        
    # displacement induced by a fault from depth b to the surface
    ux1,uy1=ds_eq_2d(x,b,dip,s)
    # displacement induced by a fault from depth t to the surface
    if t>0:
        ux2,uy2=ds_eq_2d(x,t,dip,s)
    else:
        ux2,uy2=ux1*0.0,uy1*0.0
        
    ux=ux1-ux2
    uy=uy1-uy2
    
    return(ux,uy)


###################################################################
def ds_ramp_eq_2d(x,x1,y1,dip1,l1,D,y2,s,s2):
###################################################################
    """
    Calculates horizontal and vertical displacements for a 2D ramp made of of two faults\
    fault_1 has coordinates (x1,y1) with dip dip1 and length l1
    fault_2 has coordinates starting at (x1,y1) and ending at (x1+D,y2) with y2 < y1   

    :param x     : numpy array of location of points (horizontal distance from the fault top fault edge)
    :param x1,y1,dip1,l1,D,y2: fault coordinates
    :param s     : slip on the fault (s positive is a reverse fault)
    :return ux,uy: two 1D numpy arrays of horizontal and vertical prediction displacements
    :note: dip1 controls the ramp dip orientation. dip1 > 0 means dipping along positive x-axis (to the right). dip1 = 0 will raise an error.
    """

    import numpy as np

    if y2 > y1 :
        raise ValueError(('!!! ERROR: y2 must be lesser than y1: %lf > %lf' % (y2,y1)))

    # fault 1

    top1 = y1
    bottom1 = y1 + np.sqrt( np.sin( np.radians(dip1) )**2 ) * l1
    
    ux1,uy1 = ds_buried_eq_2d(x-x1,top1,bottom1,dip1,s,ref_x='true_loc')

    # fault 2
    
    if (D > 0) and (y1 > y2) :
        
        x2      = x1 + D
        
        top2    = y2 
        bottom2 = y1 
        
        abs_dip2 = np.degrees( np.arctan( ( bottom2 - top2 ) / D ) ) 
        
        #s2      = s / np.cos(np.radians(abs_dip2))
        
        dip2 = np.sign(dip1) * abs_dip2
        
        ux2,uy2 = ds_buried_eq_2d(x-x2,top2,bottom2,dip2,s2,ref_x='true_loc')

    else:
        ux2 = ux1 * 0.0
        uy2 = uy1 * 0.0

    # ramp = fault 1 + fault 2

    ux=ux1+ux2
    uy=uy1+uy2

    return(ux,uy)

    
"""
Useful functions for two-dimensional models on a profile
"""


###################################################################
def unit_vector_profile(ilon,ilat,flon,flat):
###################################################################

    """
    Returns the unit direction vector in ENU coordinates for a profile
    
    :param ilon: longitude in dec.deg for initial point of the profile
    :param ilat: latitude in dec.deg for initial point of the profile
    :param flon: longitude in dec.deg for end point of the profile
    :param flat: latitude in dec.deg for end point of the profile
    :param N   : unit vector in ENU coordinate system as a 1D numpy array
    """
    
    # import
    
    import numpy as np
    
    from pyacs.lib import coordinates as Coordinates
    from pyacs.lib.gmtpoint import GMT_Point
    import pyacs.lib.vectors as vectors
    import pyacs.lib.euler as euler

    (xi,yi,zi)=Coordinates.geo2xyz(np.radians(ilon),np.radians(ilat),0.0)
    Xi=np.array([xi,yi,zi])
    MI=GMT_Point(lon=ilon,lat=ilat)
    
    (xs,ys,zs)=Coordinates.geo2xyz(np.radians(flon),np.radians(flat),0.0)
    Xs=np.array([xs,ys,zs])
    MS=GMT_Point(lon=flon,lat=flat)

    
    POLE=vectors.vector_product(Xi,Xs)
    POLE=vectors.normalize(POLE)
    
    (llambda,phi,omega)=euler.rot2euler(POLE[0],POLE[1],POLE[2])
    
    Mmiddle=GMT_Point(lon=(ilon+flon)/2.,lat=(ilat+flat)/2.)
    
    (x,y,z)=Coordinates.geo2xyz(np.radians(Mmiddle.lon),np.radians(Mmiddle.lat),0.0)
    
    R=Coordinates.mat_rot_general_to_local(np.radians(Mmiddle.lon),np.radians(Mmiddle.lat),Mmiddle.he)
    OM=np.array([x,y,z])
    OM=vectors.normalize(OM)
    
    unit_parallele_xyz=vectors.vector_product(POLE,OM)
    unit_parallele_enu=np.dot(R,unit_parallele_xyz)

    return(unit_parallele_enu)

###################################################################
def ds_2_los(ux,uy,u_profile,u_los):
###################################################################
    """
    Returns the LOS displacement given the unit direction vector of a profile and LOS
    
    :param ux,uy: 1D numpy array including horizontal and vertical displacements
    :param u_profile: unit direction vector in ENU coordinates
    :param u_los : unit direction vector for LOS in ENU coordinates
    :return dlos : 1D numpy array of displacement in the LOS
    """
    import numpy as np
    
    # horizontal displacement in EN
    DIS_ENU=ux.reshape(-1,1)*u_profile
    # Up displacement
    DIS_ENU[:,2]=uy
    
    LOS=np.dot(DIS_ENU,u_los)

    
    
    return(LOS)





















