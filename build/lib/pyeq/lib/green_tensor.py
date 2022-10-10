"""
Manipulate the Green tensor
"""

def GREEN_ENU_TO_GREEN_RAKE_MAIN_RAKE_CONJUGATE(GREEN, GEOMETRY, SGEOMETRY, rake_type, rake_or_pole):
    """
    GREEN IS A TENSOR OF DIM 4
    The input GREEN TENSOR is
    GREEN(i,j,k,l) is the prediction for dislocation i at site j component k for rake l
    k=0,1,2 = east, north, up
    l=0,1 : rake_00 & rake_90
    I want for this code
    NEW_GREEN(i,j,k,l) is the prediction for dislocation j at site i component k for rake l
    k=0,1,2 = east, north, up
    l=0,1 : rake_principal & rake_conjugate
    
    the rake can be defined either by an euler pole and geometry (providing position strike and dip to calculate the relative slip) or a fixed rake
    
    """

    import numpy as np
    import pyacs.lib.faultslip
    import pyacs.lib.euler
     
    # dealing with the principal rake
    
    RAKE=np.zeros(GEOMETRY.shape[0])
    VEL_FROM_EULER=np.zeros(GEOMETRY.shape[0])
    
    if rake_type.lower() =='euler':
        
        _tmp,elon,elat,ew,_motion_type=rake_or_pole.split('/')
#        pole=("%s/%s/%s" % (elon,elat,ew))
    
        for i in range(GEOMETRY.shape[0]):
            [x,y,strike,dip]=[SGEOMETRY.centroid_long[i],SGEOMETRY.centroid_lat[i],SGEOMETRY.strike[i],SGEOMETRY.dip[i]]
    
            RAKE[i]=pyacs.lib.faultslip.rake_from_euler(x,y,strike, dip, rake_or_pole)
            (ve,vn) = pyacs.lib.euler.vel_from_euler( x , y , float(elon) , float(elat) , float(ew) )
            VEL_FROM_EULER[i] = np.sqrt( ve**2 + vn**2 )
    else:
        RAKE=RAKE+float(rake_or_pole)

    print('gt RAKE')
    print(RAKE.shape)



    RAKE_RADIANS=np.radians(RAKE)
    CONJUGATE_RAKE_RADIANS=np.radians(RAKE+90.0)
    
    
    
    NEW_GREEN=np.zeros((GREEN.shape[1],GREEN.shape[0],GREEN.shape[2],GREEN.shape[3]))
    
    GREEN_4GPS_EAST_RAKE_00=GREEN[:,:,0,0].T
    GREEN_4GPS_EAST_RAKE_90=GREEN[:,:,0,1].T
    GREEN_4GPS_NORTH_RAKE_00=GREEN[:,:,1,0].T
    GREEN_4GPS_NORTH_RAKE_90=GREEN[:,:,1,1].T
    GREEN_4UP_RAKE_00=GREEN[:,:,2,0].T
    GREEN_4UP_RAKE_90=GREEN[:,:,2,1].T

    print('GREEN_4GPS_EAST_RAKE_00')
    print(GREEN_4GPS_EAST_RAKE_00.shape)

    # Now calculating the Green's functions in the principal rake direction
    
    GREEN_4GPS_EAST_RAKE_PRINCIPAL =np.cos(RAKE_RADIANS)*GREEN_4GPS_EAST_RAKE_00+np.sin(RAKE_RADIANS)*GREEN_4GPS_EAST_RAKE_90
    GREEN_4GPS_NORTH_RAKE_PRINCIPAL=np.cos(RAKE_RADIANS)*GREEN_4GPS_NORTH_RAKE_00+np.sin(RAKE_RADIANS)*GREEN_4GPS_NORTH_RAKE_90
    
    GREEN_4GPS_EAST_RAKE_CONJUGATE=np.cos(CONJUGATE_RAKE_RADIANS)*GREEN_4GPS_EAST_RAKE_00+np.sin(CONJUGATE_RAKE_RADIANS)*GREEN_4GPS_EAST_RAKE_90
    GREEN_4GPS_NORTH_RAKE_CONJUGATE=np.cos(CONJUGATE_RAKE_RADIANS)*GREEN_4GPS_NORTH_RAKE_00+np.sin(CONJUGATE_RAKE_RADIANS)*GREEN_4GPS_NORTH_RAKE_90
    
    GREEN_4GPS_UP_RAKE_PRINCIPAL=np.cos(RAKE_RADIANS)*GREEN_4UP_RAKE_00+np.sin(RAKE_RADIANS)*GREEN_4UP_RAKE_90
    GREEN_4GPS_UP_RAKE_CONJUGATE=np.cos(CONJUGATE_RAKE_RADIANS)*GREEN_4UP_RAKE_00+np.sin(CONJUGATE_RAKE_RADIANS)*GREEN_4UP_RAKE_90
    
    NEW_GREEN[:,:,0,0]=GREEN_4GPS_EAST_RAKE_PRINCIPAL
    NEW_GREEN[:,:,1,0]=GREEN_4GPS_NORTH_RAKE_PRINCIPAL
    NEW_GREEN[:,:,2,0]=GREEN_4GPS_UP_RAKE_PRINCIPAL
    
    NEW_GREEN[:,:,0,1]=GREEN_4GPS_EAST_RAKE_CONJUGATE
    NEW_GREEN[:,:,1,1]=GREEN_4GPS_NORTH_RAKE_CONJUGATE
    NEW_GREEN[:,:,2,1]=GREEN_4GPS_UP_RAKE_CONJUGATE

    return(NEW_GREEN)