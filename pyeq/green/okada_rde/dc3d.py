"""
This a a pure python direct translation from Okada's dc3d fortran original routine.
ref: Okada (1992) [Bull. Seism. Soc. Am., 82, 1018-1040]
http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html
http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html
http://www.bosai.go.jp/study/application/dc3d/download/DC3Dfortran.txt
"""

# global variables

f0 = 0.
f1 = 1.
f2 = 2.
f3 = 3.
f4 = 4.
f5 = 5.
f7 = 7.
f8 = 8.
f9 = 9.
f10 = 10.
f15 = 15.

pi2 = 6.283185307179586                                     
eps = 1.E-6


###############################################################################
def dc3d(alpha, \
###############################################################################
         x, y, z, \
         depth, \
         dip, \
         al1, al2, \
         aw1, aw2, \
         disl1, disl2, disl3):
    """                                                                       
    c********************************************************************   
    c*****                                                          *****   
    c*****    displacement and strain at depth                      *****   
    c*****    due to buried finite fault in a semiinfinite medium   *****   
    c*****              coded by  y.okada_rde ... sep.1991              *****
    c*****              revised ... nov.1991, apr.1992, may.1993,   *****   
    c*****                          jul.1993                        *****   
    c********************************************************************   
    c                                                                       
    c***** input                                                            
    c*****   alpha : medium constant  (lambda+myu)/(lambda+2*myu)           
    c*****   x,y,z : coordinate of observing point                          
    c*****   depth : depth of reference point                               
    c*****   dip   : dip-angle (degree)                                     
    c*****   al1,al2   : fault length range                                 
    c*****   aw1,aw2   : fault width range                                  
    c*****   disl1-disl3 : strike-, dip-, tensile-dislocations              
    c                                                                       
    c***** output                                                           
    c*****   ux, uy, uz  : displacement ( unit=(unit of disl)               
    c*****   uxx,uyx,uzx : x-derivative ( unit=(unit of disl) /             
    c*****   uxy,uyy,uzy : y-derivative        (unit of x,y,z,depth,al,aw) )
    c*****   uxz,uyz,uzz : z-derivative                                     
    c*****   iret        : return code  ( =0....normal,   =1....singular )  
    """

#    global p, q, s, t, xy, x2, y2, d2, r, r2, r3, r5, qr, qrx, a3, a5, b3, c3, uy, vy, wy, uz, vz, wz
    import numpy as np
    
    u = np.zeros(12) + f0
    du = np.zeros(12) + f0
    dua = np.zeros(12) + f0
    dub = np.zeros(12) + f0
    duc = np.zeros(12) + f0
    
    xi = np.zeros(2)                                                            
    et = np.zeros(2)                                                            
    kxi = np.zeros(2)                                                            
    ket = np.zeros(2)                                                            
#      common /c0/dummy(5),sd,cd,dumm(5)                                 
#      dimension  xi(2),et(2),kxi(2),ket(2)                              
#      dimension  u(12),du[12),dua(12),dub(12),duc(12)                   
#      data  f0,eps/ 0.d0, 1.d-6 /                                       

    if (z  > 0.):
        print(' ** positive z was given in sub-dc3d') 
        print(z)
        sys.exit()
    aalpha = alpha                                                      
    ddip = dip                                                          

    dccon0(aalpha, ddip)                                          

    zz = z                                                              
    dd1 = disl1                                                         
    dd2 = disl2                                                         
    dd3 = disl3                                                         
    xi[0] = x - al1                                                       
    xi[1] = x - al2                                                       
    if (np.abs(xi[0]) <= eps): xi[0] = f0                                   
    if (np.abs(xi[1]) <= eps): xi[1] = f0
#======================================                                 
#=====  real-source contribution  =====                                 
#======================================                                 
    
    d = depth + z                                                         
    p = y * cd + d * sd                                                       
    q = y * sd - d * cd    
    
                                                       
    et[0] = p - aw1                                                       
    et[1] = p - aw2                                                       
    
    if (np.abs(q) < eps):  q = f0                                          
    if (np.abs(et[0]) < eps): et[0] = f0                                   
    if (np.abs(et[1]) < eps): et[1] = f0                                   

#--------------------------------                                       
#----- reject singular case -----                                       
#--------------------------------                                       
#----- on fault edge                                                    
    if ((q == f0) and ((xi[0] * xi[1] <= f0 and et[0] * et[1] == f0) or (et[0] * et[1] <= f0 and xi[0] * xi[1] == f0))):
        return np.zeros(12)                                                       
#----- on negative extension of fault edge                              
    kxi[0] = 0                                                          
    kxi[1] = 0                                                          
    ket[0] = 0                                                          
    ket[1] = 0                                                          
    r12 = np.sqrt(xi[0] * xi[0] + et[1] * et[1] + q * q)                            
    r21 = np.sqrt(xi[1] * xi[1] + et[0] * et[0] + q * q)                            
    r22 = np.sqrt(xi[1] * xi[1] + et[1] * et[1] + q * q)                            
    if (xi[0] < f0 and r21 + xi[1] < eps): kxi[0] = 1                   
    if (xi[0] < f0 and r22 + xi[1] < eps): kxi[1] = 1                   
    if (et[0] < f0 and r12 + et[1] < eps): ket[0] = 1                   
    if (et[0] < f0 and r22 + et[1] < eps): ket[1] = 1
#=====                                                                  

    for k in [0, 1]:                                                      
        for j in [0, 1]:                                                      
            dccon2(xi[j], et[k], q, sd, cd, kxi[k], ket[j])                  
            dua = ua(xi[j], et[k], q, dd1, dd2, dd3)                          
#-----                                                                  
            for i in [0, 3, 6, 9]:
                du[i] = -dua[i]                                               
                du[i + 1] = -dua[i + 1] * cd + dua[i + 2] * sd                              
                du[i + 2] = -dua[i + 1] * sd - dua[i + 2] * cd                              

            du[ 9] = -du[9]                                                
            du[10] = -du[10]                                              
            du[11] = -du[11]                                              

            for i in np.arange(12):
                if ((j + k) != 1): u[i] = u[i] + du[i]                                  
                if ((j + k) == 1): u[i] = u[i] - du[i]                                  
    
#=======================================                                
#=====  image-source contribution  =====                                
#=======================================                                
    d = depth - z                                                         
    p = y * cd + d * sd                                                       
    q = y * sd - d * cd                                                       
    et[0] = p - aw1                                                       
    et[1] = p - aw2                                                       
    if (np.abs(q) < eps):  q = f0                                          
    if (np.abs(et[0]) < eps): et[0] = f0                                   
    if (np.abs(et[1]) < eps): et[1] = f0
#--------------------------------                                       
#----- reject singular case -----                                       
#--------------------------------                                       
#----- on fault edge                                                    
    if ((q == f0) and ((xi[0] * xi[1] <= f0 and et[0] * et[1] == f0) or (et[0] * et[1] <= f0 and xi[0] * xi[1] == f0))):       
        return np.zeros(12)                                                       
#----- on negative extension of fault edge                              
    kxi[0] = 0                                                          
    kxi[1] = 0                                                          
    ket[0] = 0                                                          
    ket[1] = 0                                                          
    r12 = np.sqrt(xi[0] * xi[0] + et[1] * et[1] + q * q)                            
    r21 = np.sqrt(xi[1] * xi[1] + et[0] * et[0] + q * q)                            
    r22 = np.sqrt(xi[1] * xi[1] + et[1] * et[1] + q * q)                            
    if (xi[0] < f0 and r21 + xi[1] < eps): kxi[0] = 1                   
    if (xi[0] < f0 and r22 + xi[1] < eps): kxi[1] = 1                   
    if (et[0] < f0 and r12 + et[1] < eps): ket[0] = 1                   
    if (et[0] < f0 and r22 + et[1] < eps): ket[1] = 1
#=====                                                                  
    for k in [0, 1]:  # do 334 k=1,2                                                      
        for j in [0, 1]:  # do 333 j=1,2                                                      
            dccon2(xi[j], et[k], q, sd, cd, kxi[k], ket[j])                  
            dua = ua(xi[j], et[k], q, dd1, dd2, dd3)                          
            dub = ub(xi[j], et[k], q, dd1, dd2, dd3)                          
            duc = uc(xi[j], et[k], q, zz, dd1, dd2, dd3)                       
#-----                                                                  
            for i in [0, 3, 6, 9]:  # do 330 i=1,10,3                                                 
                du[i] = dua[i] + dub[i] + z * duc[i]                                  
                du[i + 1] = (dua[i + 1] + dub[i + 1] + z * duc[i + 1]) * cd - (dua[i + 2] + dub[i + 2] + z * duc[i + 2]) * sd                     
                du[i + 2] = (dua[i + 1] + dub[i + 1] - z * duc[i + 1]) * sd + (dua[i + 2] + dub[i + 2] - z * duc[i + 2]) * cd                     

            du[9] = du[9] + duc[0]                                          
            du[10] = du[10] + duc[1] * cd - duc[2] * sd                             
            du[11] = du[11] - duc[1] * sd - duc[2] * cd                             

            for i in np.arange(12):
                if (j + k != 1): u[i] = u[i] + du[i]                                  
                if (j + k == 1): u[i] = u[i] - du[i]                                  


    return u


###############################################################################
def ua(xi, et, q, disl1, disl2, disl3):                       
###############################################################################
    """"
    c********************************************************************   
    c*****    displacement and strain at depth (part-a)             *****   
    c*****    due to buried finite fault in a semiinfinite medium   *****   
    c********************************************************************   
    c                                                                       
    c***** input                                                            
    c*****   xi,et,q : station coordinates in fault system                  
    c*****   disl1-disl3 : strike-, dip-, tensile-dislocations              
    c***** output                                                           
    c*****   u(12) : displacement and their derivatives                     
    """

    import numpy as np

    du = np.zeros(12)
    u = np.zeros(12)
                                                                    
#      common /c0/alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d  
#      common /c2/xi2,et2,q2,r,r2,r3,r5,y,d,tt,alx,ale,x11,y11,x32,y32,  
#     *           ey,ez,fy,fz,gy,gz,hy,hz                                
#      data f0,f2,pi2/0.d0,2.d0,6.283185307179586d0/                     
#-----                                                                  
    xy = xi * y11                                                         
    qx = q * x11                                                         
    qy = q * y11                                                         

#======================================                                 
#=====  strike-slip contribution  =====                                 
#======================================                                 
    if (disl1 != f0):                                              
        du[ 0] = tt / f2 + alp2 * xi * qy                                    
        du[ 1] = alp2 * q / r                                      
        du[ 2] = alp1 * ale - alp2 * q * qy                                     
        du[ 3] = -alp1 * qy - alp2 * xi2 * q * y32                                
        du[ 4] = -alp2 * xi * q / r3                                  
        du[ 5] = alp1 * xy + alp2 * xi * q2 * y32                                
        du[ 6] = alp1 * xy * sd + alp2 * xi * fy + d / f2 * x11                  
        du[ 7] = alp2 * ey                              
        du[ 8] = alp1 * (cd / r + qy * sd) - alp2 * q * fy                            
        du[ 9] = alp1 * xy * cd + alp2 * xi * fz + y / f2 * x11                  
        du[10] = alp2 * ez                              
        du[11] = -alp1 * (sd / r - qy * cd) - alp2 * q * fz                            
        
        u = u + disl1 / pi2 * du

#======================================                                 
#=====    dip-slip contribution   =====                                 
#======================================                                 
    if (disl2 != f0):                                              
        du[ 0] = alp2 * q / r                                      
        du[ 1] = tt / f2 + alp2 * et * qx                                    
        du[ 2] = alp1 * alx - alp2 * q * qx                                     
        du[ 3] = -alp2 * xi * q / r3                                    
        du[ 4] = -qy / f2 - alp2 * et * q / r3                                    
        du[ 5] = alp1 / r + alp2 * q2 / r3                                      
        du[ 6] = alp2 * ey                            
        du[ 7] = alp1 * d * x11 + xy / f2 * sd + alp2 * et * gy                         
        du[ 8] = alp1 * y * x11 - alp2 * q * gy                          
        du[ 9] = alp2 * ez                            
        du[10] = alp1 * y * x11 + xy / f2 * cd + alp2 * et * gz                         
        du[11] = -alp1 * d * x11 - alp2 * q * gz                          

        u = u + disl2 / pi2 * du

#========================================                               
#=====  tensile-fault contribution  =====                               
#========================================                               
    if (disl3 != f0):                                              
        du[ 0] = -alp1 * ale - alp2 * q * qy                                     
        du[ 1] = -alp1 * alx - alp2 * q * qx                                     
        du[ 2] = tt / f2 - alp2 * (et * qx + xi * qy)                            
        du[ 3] = -alp1 * xy + alp2 * xi * q2 * y32                                
        du[ 4] = -alp1 / r + alp2 * q2 / r3                                    
        du[ 5] = -alp1 * qy - alp2 * q * q2 * y32                                 
        du[ 6] = -alp1 * (cd / r + qy * sd) - alp2 * q * fy                           
        du[ 7] = -alp1 * y * x11 - alp2 * q * gy                           
        du[ 8] = alp1 * (d * x11 + xy * sd) + alp2 * q * hy                           
        du[ 9] = alp1 * (sd / r - qy * cd) - alp2 * q * fz                           
        du[10] = alp1 * d * x11 - alp2 * q * gz                           
        du[11] = alp1 * (y * x11 + xy * cd) + alp2 * q * hz                           

        u = u + disl3 / pi2 * du

    return u


###############################################################################
def ub(xi, et, q, disl1, disl2, disl3):                    
###############################################################################
    """                                                                       
    c********************************************************************   
    c*****    displacement and strain at depth (part-b)             *****   
    c*****    due to buried finite fault in a semiinfinite medium   *****   
    c********************************************************************   
    c                                                                       
    c***** input                                                            
    c*****   xi,et,q : station coordinates in fault system                  
    c*****   disl1-disl3 : strike-, dip-, tensile-dislocations              
    c***** output                                                           
    c*****   u(12) : displacement and their derivatives                     
    """

    import numpy as np
      
    du = np.zeros(12)  
                                                          
#      common /c0/alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d  
#      common /c2/xi2,et2,q2,r,r2,r3,r5,y,d,tt,alx,ale,x11,y11,x32,y32,  
#     *           ey,ez,fy,fz,gy,gz,hy,hz                                
#      data  f0,f1,f2,pi2/0.d0,1.d0,2.d0,6.283185307179586d0/            
#-----
    rd = r + d                                                            
    d11 = f1 / (r * rd)                                                     
    aj2 = xi * y / rd * d11                                                   
    aj5 = -(d + y * y / rd) * d11                                               

    if (cd != f0):
        if (xi == f0):
            ai4 = f0                                                        
        else:                                                            
            x = np.sqrt(xi2 + q2)                                               
            ai4 = f1 / cdcd * (xi / rd * sdcd + f2 * np.arctan((et * (x + q * cd) + x * (r + x) * sd) / (xi * (r + x) * cd)))        
                                                                 
        ai3 = (y * cd / rd - ale + sd * np.log(rd)) / cdcd                              
        ak1 = xi * (d11 - y11 * sd) / cd                                          
        ak3 = (q * y11 - y * d11) / cd                                            
        aj3 = (ak1 - aj2 * sd) / cd                                             
        aj6 = (ak3 - aj5 * sd) / cd                                             
    else:
        rd2 = rd * rd                                                       
        ai3 = (et / rd + y * q / rd2 - ale) / f2                                      
        ai4 = xi * y / rd2 / f2                                                 
        ak1 = xi * q / rd * d11                                                 
        ak3 = sd / rd * (xi2 * d11 - f1)                                          
        aj3 = -xi / rd2 * (q2 * d11 - f1 / f2)                                      
        aj6 = -y / rd2 * (xi2 * d11 - f1 / f2)                                      

#-----                                                                  
    xy = xi * y11                                                         
    ai1 = -xi / rd * cd - ai4 * sd                                              
    ai2 = np.log(rd) + ai3 * sd                                              
    ak2 = f1 / r + ak3 * sd                                                  
    ak4 = xy * cd - ak1 * sd                                                 
    aj1 = aj5 * cd - aj6 * sd                                                
    aj4 = -xy - aj2 * cd + aj3 * sd                                             
#=====                                                                  
    u = np.zeros(12)
    qx = q * x11                                                          
    qy = q * y11                                                          

#======================================                                 
#=====  strike-slip contribution  =====                                 
#======================================                                 
    if (disl1 != f0):                                              
        du[ 0] = -xi * qy - tt - alp3 * ai1 * sd                                   
        du[ 1] = -q / r + alp3 * y / rd * sd                                  
        du[ 2] = q * qy - alp3 * ai2 * sd                                   
        du[ 3] = xi2 * q * y32 - alp3 * aj1 * sd                                  
        du[ 4] = xi * q / r3 - alp3 * aj2 * sd                                  
        du[ 5] = -xi * q2 * y32 - alp3 * aj3 * sd                                  
        du[ 6] = -xi * fy - d * x11 + alp3 * (xy + aj4) * sd                           
        du[ 7] = -ey + alp3 * (f1 / r + aj5) * sd                         
        du[ 8] = q * fy - alp3 * (qy - aj6) * sd                           
        du[ 9] = -xi * fz - y * x11 + alp3 * ak1 * sd                                
        du[10] = -ez + alp3 * y * d11 * sd                              
        du[11] = q * fz + alp3 * ak2 * sd                                

        u = u + disl1 / pi2 * du                                       

#======================================                                 
#=====    dip-slip contribution   =====                                 
#======================================                                 
    if (disl2 != f0):
        du[ 0] = -q / r + alp3 * ai3 * sdcd                                 
        du[ 1] = -et * qx - tt - alp3 * xi / rd * sdcd                               
        du[ 2] = q * qx + alp3 * ai4 * sdcd                                 
        du[ 3] = xi * q / r3 + alp3 * aj4 * sdcd                              
        du[ 4] = et * q / r3 + qy + alp3 * aj5 * sdcd                              
        du[ 5] = -q2 / r3 + alp3 * aj6 * sdcd                              
        du[ 6] = -ey + alp3 * aj1 * sdcd                              
        du[ 7] = -et * gy - xy * sd + alp3 * aj2 * sdcd                              
        du[ 8] = q * gy + alp3 * aj3 * sdcd                              
        du[ 9] = -ez - alp3 * ak3 * sdcd                              
        du[10] = -et * gz - xy * cd - alp3 * xi * d11 * sdcd                           
        du[11] = q * gz - alp3 * ak4 * sdcd                              

        u = u + disl2 / pi2 * du                                       
#========================================                               
#=====  tensile-fault contribution  =====                               
#========================================                               
    if (disl3 != f0):
        du[ 0] = q * qy - alp3 * ai3 * sdsd                           
        du[ 1] = q * qx + alp3 * xi / rd * sdsd                         
        du[ 2] = et * qx + xi * qy - tt - alp3 * ai4 * sdsd                           
        du[ 3] = -xi * q2 * y32 - alp3 * aj4 * sdsd                                
        du[ 4] = -q2 / r3 - alp3 * aj5 * sdsd                                
        du[ 5] = q * q2 * y32 - alp3 * aj6 * sdsd                                
        du[ 6] = q * fy - alp3 * aj1 * sdsd                                     
        du[ 7] = q * gy - alp3 * aj2 * sdsd                                     
        du[ 8] = -q * hy - alp3 * aj3 * sdsd                                     
        du[ 9] = q * fz + alp3 * ak3 * sdsd                                     
        du[10] = q * gz + alp3 * xi * d11 * sdsd                                  
        du[11] = -q * hz + alp3 * ak4 * sdsd                                     

        u = u + disl3 / pi2 * du                                       
      
    return u                                                            


###############################################################################
def uc(xi, et, q, z, disl1, disl2, disl3):                     
###############################################################################
    """                                                                       
    c********************************************************************   
    c*****    displacement and strain at depth (part-c)             *****   
    c*****    due to buried finite fault in a semiinfinite medium   *****   
    c********************************************************************   
    c                                                                       
    c***** input                                                            
    c*****   xi,et,q,z   : station coordinates in fault system              
    c*****   disl1-disl3 : strike-, dip-, tensile-dislocations              
    c***** output                                                           
    c*****   u(12) : displacement and their derivatives                     
    """
    
    import numpy as np
                                          
    du = np.zeros(12)
    
#      common /c0/alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d  
#      common /c2/xi2,et2,q2,r,r2,r3,r5,y,d,tt,alx,ale,x11,y11,x32,y32,  
#     *           ey,ez,fy,fz,gy,gz,hy,hz                                
#      data f0,f1,f2,f3,pi2/0.d0,1.d0,2.d0,3.d0,6.283185307179586d0/     
#-----                                                                  
    c = d + z                                                             
    x53 = (8.*r2 + 9.*r * xi + f3 * xi2) * x11 * x11 * x11 / r2                     
    y53 = (8.*r2 + 9.*r * et + f3 * et2) * y11 * y11 * y11 / r2                     
    h = q * cd - z                                                          
    z32 = sd / r3 - h * y32                                                   
    z53 = f3 * sd / r5 - h * y53                                                
    y0 = y11 - xi2 * y32                                                    
    z0 = z32 - xi2 * z53                                                    
    ppy = cd / r3 + q * y32 * sd                                                
    ppz = sd / r3 - q * y32 * cd                                                
    qq = z * y32 + z32 + z0                                                   
    qqy = f3 * c * d / r5 - qq * sd                                               
    qqz = f3 * c * y / r5 - qq * cd + q * y32                                         
    xy = xi * y11                                                         
    qx = q * x11                                                          
    qy = q * y11                                                          
    qr = f3 * q / r5                                                        
    cqx = c * q * x53                                                       
    cdr = (c + d) / r3                                                      
    yy0 = y / r3 - y0 * cd                                                    
#=====                                                                  
    u = np.zeros(12)

#======================================                                 
#=====  strike-slip contribution  =====                                 
#======================================                                 
    if (disl1 != f0):                                              
        du[ 0] = alp4 * xy * cd - alp5 * xi * q * z32                     
        du[ 1] = alp4 * (cd / r + f2 * qy * sd) - alp5 * c * q / r3                       
        du[ 2] = alp4 * qy * cd - alp5 * (c * et / r3 - z * y11 + xi2 * z32)      
        du[ 3] = alp4 * y0 * cd - alp5 * q * z0                  
        du[ 4] = -alp4 * xi * (cd / r3 + f2 * q * y32 * sd) + alp5 * c * xi * qr               
        du[ 5] = -alp4 * xi * q * y32 * cd + alp5 * xi * (f3 * c * et / r5 - qq)    
        du[ 6] = -alp4 * xi * ppy * cd - alp5 * xi * qqy                          
        du[ 7] = alp4 * f2 * (d / r3 - y0 * sd) * sd - y / r3 * cd - alp5 * (cdr * sd - et / r3 - c * y * qr)           
        du[ 8] = -alp4 * q / r3 + yy0 * sd + alp5 * (cdr * cd + c * d * qr - (y0 * cd + q * z0) * sd) 
        du[ 9] = alp4 * xi * ppz * cd - alp5 * xi * qqz                          
        du[10] = alp4 * f2 * (y / r3 - y0 * cd) * sd + d / r3 * cd - alp5 * (cdr * cd + c * d * qr)   
        du[11] = yy0 * cd - alp5 * (cdr * sd - c * y * qr - y0 * sdsd + q * z0 * cd) 
        
        u = u + disl1 / pi2 * du                                       
#======================================                                 
#=====    dip-slip contribution   =====                                 
#======================================                                 
    if (disl2 != f0):                                              
        du[ 0] = alp4 * cd / r - qy * sd - alp5 * c * q / r3                           
        du[ 1] = alp4 * y * x11 - alp5 * c * et * q * x32                       
        du[ 2] = -d * x11 - xy * sd - alp5 * c * (x11 - q2 * x32)                   
        du[ 3] = -alp4 * xi / r3 * cd + alp5 * c * xi * qr + xi * q * y32 * sd                
        du[ 4] = -alp4 * y / r3 + alp5 * c * et * qr                             
        du[ 5] = d / r3 - y0 * sd + alp5 * c / r3 * (f1 - f3 * q2 / r2)                  
        du[ 6] = -alp4 * et / r3 + y0 * sdsd - alp5 * (cdr * sd - c * y * qr)                
        du[ 7] = alp4 * (x11 - y * y * x32) - alp5 * c * ((d + f2 * q * cd) * x32 - y * et * q * x53) 
        du[ 8] = xi * ppy * sd + y * d * x32 + alp5 * c * ((y + f2 * q * sd) * x32 - y * q2 * x53)   
        du[ 9] = -q / r3 + y0 * sdcd - alp5 * (cdr * cd + c * d * qr)                
        du[10] = alp4 * y * d * x32 - alp5 * c * ((y - f2 * q * sd) * x32 + d * et * q * x53) 
        du[11] = -xi * ppz * sd + x11 - d * d * x32 - alp5 * c * ((d - f2 * q * cd) * x32 - d * q2 * x53) 

        u = u + disl2 / pi2 * du
                                               
#========================================                               
#=====  tensile-fault contribution  =====                               
#========================================                               
    if (disl3 != f0):                                              
        du[ 0] = -alp4 * (sd / r + qy * cd) - alp5 * (z * y11 - q2 * z32)                
        du[ 1] = alp4 * f2 * xy * sd + d * x11 - alp5 * c * (x11 - q2 * x32)                
        du[ 2] = alp4 * (y * x11 + xy * cd) + alp5 * q * (c * et * x32 + xi * z32)           
        du[ 3] = alp4 * xi / r3 * sd + xi * q * y32 * cd + alp5 * xi * (f3 * c * et / r5 - f2 * z32 - z0)
        du[ 4] = alp4 * f2 * y0 * sd - d / r3 + alp5 * c / r3 * (f1 - f3 * q2 / r2)             
        du[ 5] = -alp4 * yy0 - alp5 * (c * et * qr - q * z0)                 
        du[ 6] = alp4 * (q / r3 + y0 * sdcd) + alp5 * (z / r3 * cd + c * d * qr - q * z0 * sd)    
        du[ 7] = -alp4 * f2 * xi * ppy * sd - y * d * x32 + alp5 * c * ((y + f2 * q * sd) * x32 - y * q2 * x53)            
        du[ 8] = -alp4 * (xi * ppy * cd - x11 + y * y * x32) + alp5 * (c * ((d + f2 * q * cd) * x32 - y * et * q * x53) + xi * qqy) 
        du[ 9] = -et / r3 + y0 * cdcd - alp5 * (z / r3 * sd - c * y * qr - y0 * sdsd + q * z0 * cd)  
        du[10] = alp4 * f2 * xi * ppz * sd - x11 + d * d * x32 - alp5 * c * ((d - f2 * q * cd) * x32 - d * q2 * x53)            
        du[11] = alp4 * (xi * ppz * cd + y * d * x32) + alp5 * (c * ((y - f2 * q * sd) * x32 + d * et * q * x53) + xi * qqz) 

        u = u + disl3 / pi2 * du
    
    return u


###############################################################################
def  dccon0(alpha, dip):                                     
###############################################################################
    """                                                                     
    c*******************************************************************    
    c*****   calculate medium constants and fault-dip constants    *****    
    c*******************************************************************    
    c                                                                       
    c***** input                                                            
    c*****   alpha : medium constant  (lambda+myu)/(lambda+2*myu)           
    c*****   dip   : dip-angle (degree)                                     
    c### caution ### if cos(dip) is sufficiently small, it is set to zero   
    """                                                                     
    import numpy as np
    global alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d, c2d
    
#      common /c0/alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d  
#      data f0,f1,f2,pi2/0.d0,1.d0,2.d0,6.283185307179586d0/             
#      data eps/1.d-6/                                                   

    alp1 = (f1 - alpha) / f2                                                
    alp2 = alpha / f2                                                    
    alp3 = (f1 - alpha) / alpha                                             
    alp4 = f1 - alpha                                                    
    alp5 = alpha                                                       

    p18 = pi2 / 360.                                                    
    sd = np.sin(dip * p18)                                                  
    cd = np.cos(dip * p18)                                                  

    if (np.abs(cd) <= eps) :                                          
        cd = f0                                                           
        if (sd >= f0): sd = f1                                             
        if (sd <= f0): sd = -f1                                             

    sdsd = sd * sd                                                        
    cdcd = cd * cd                                                        
    sdcd = sd * cd                                                        
    s2d = f2 * sdcd                                                       
    c2d = cdcd - sdsd                                                     
    
    return 


###############################################################################
def dccon1(x, y, d):
###############################################################################
    """
    c********************************************************************** 
    c*****   calculate station geometry constants for point source    ***** 
    c********************************************************************** 
    c                                                                       
    c***** input                                                            
    c*****   x,y,d : station coordinates in fault system                    
    c### caution ### if x,y,d are sufficiently small, they are set to zero  
    
    """
    
    import numpy as np
    global p, q, s, t, xy, x2, y2, d2, r, r2, r3, r5, qr, qrx, a3, a5, b3, c3, uy, vy, wy, uz, vz, wz
                                                                     
#      common /c0/dummy(5),sd,cd,dumm(5)                                 
#      common /c1/p,q,s,t,xy,x2,y2,d2,r,r2,r3,r5,qr,qrx,a3,a5,b3,c3,     
#     *           uy,vy,wy,uz,vz,wz                                      
#      data  f0,f1,f3,f5,eps/0.d0,1.d0,3.d0,5.d0,1.d-6/                  

    if (np.abs(x) < eps): x = f0                                           
    if (np.abs(y) < eps): y = f0                                           
    if (np.abs(d) < eps): d = f0                                           

    p = y * cd + d * sd                                                       
    q = y * sd - d * cd                                                       
    s = p * sd + q * cd                                                       
    t = p * cd - q * sd                                                       

    xy = x * y                                                            
    x2 = x * x                                                            
    y2 = y * y                                                            
    d2 = d * d                                                            
    r2 = x2 + y2 + d2                                                       
    r = np.sqrt(r2)                                                      

    if (r == f0): 
        return                                             

    r3 = r * r2                                                          
    r5 = r3 * r2                                                          
    _r7 = r5 * r2                                                          

    a3 = f1 - f3 * x2 / r2                                                    
    a5 = f1 - f5 * x2 / r2                                                    
    b3 = f1 - f3 * y2 / r2                                                    
    c3 = f1 - f3 * d2 / r2                                                    
    
    qr = f3 * q / r5                                                        
    qrx = f5 * qr * x / r2                                                    
    
    uy = sd - f5 * y * q / r2                                                   
    uz = cd + f5 * d * q / r2                                                   
    vy = s - f5 * y * p * q / r2                                                 
    vz = t + f5 * d * p * q / r2                                                 
    wy = uy + sd                                                          
    wz = uz + cd                                                          

    return                                             


###############################################################################
def dccon2(xi, et, q, sd, cd, kxi, ket):                   
###############################################################################

    """                                                                       
    c********************************************************************** 
    c*****   calculate station geometry constants for finite source   ***** 
    c********************************************************************** 
    c                                                                       
    c***** input                                                            
    c*****   xi,et,q : station coordinates in fault system                  
    c*****   sd,cd   : sin, cos of dip-angle                                
    c*****   kxi,ket : kxi=1, ket=1 means r+xi<eps, r+et<eps, respectively  
    c                                                                       
    c### caution ### if xi,et,q are sufficiently small, they are set to zer0
    """
    
    import numpy as np
                                                                     
    global xi2, et2, q2, r, r2, r3, r5, y, d, tt, alx, ale, x11, y11, x32, y32 
    global ey, ez, fy, fz, gy, gz, hy, hz                                

#-----                                                                  
    if (np.abs(xi) < eps): xi = f0                                         
    if (np.abs(et) < eps): et = f0                                         
    if (np.abs(q) < eps):  q = f0                                         
    xi2 = xi * xi                                                         
    et2 = et * et                                                         
    q2 = q * q                                                            
    r2 = xi2 + et2 + q2                                                     
    r = np.sqrt(r2)                                                      
    if (r == f0): return                                                
    r3 = r * r2                                                          
    r5 = r3 * r2                                                          
    y = et * cd + q * sd                                                     
    d = et * sd - q * cd                                                     
#-----                                                                  
    if (q == f0): tt = f0                                                           
    else:  tt = np.arctan(xi * et / (q * r))                                           

#-----                                                                  
    if (kxi == 1):                                                 
        alx = -np.log(r - xi)                                                 
        x11 = f0                                                          
        x32 = f0                                                          
    else:                                                              
        rxi = r + xi                                                        
        alx = np.log(rxi)                                                   
        x11 = f1 / (r * rxi)                                                  
        x32 = (r + rxi) * x11 * x11 / r                                           

#-----                                                                  
    if (ket == 1):                                                 
        ale = -np.log(r - et)                                                 
        y11 = f0                                                          
        y32 = f0                                                          
    else:                                                              
        ret = r + et                                                        
        ale = np.log(ret)                                                   
        y11 = f1 / (r * ret)                                                  
        y32 = (r + ret) * y11 * y11 / r                                           

#-----                                                                  
    ey = sd / r - y * q / r3                                                    
    ez = cd / r + d * q / r3                                                    
    fy = d / r3 + xi2 * y32 * sd                                                
    fz = y / r3 + xi2 * y32 * cd                                                
    gy = f2 * x11 * sd - y * q * x32                                              
    gz = f2 * x11 * cd + d * q * x32                                              
    hy = d * q * x32 + xi * q * y32 * sd                                            
    hz = y * q * x32 + xi * q * y32 * cd                                            
