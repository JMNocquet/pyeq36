"""
This a a pure python direct translation from Okada's dc3d0 fortran original routine.
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
def ua0(x, y, d, pot1, pot2, pot3, pot4):
###############################################################################
#      common /c0/alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d  
#      common /c1/p,q,s,t,xy,x2,y2,d2,r,r2,r3,r5,qr,qrx,a3,a5,b3,c3,     
#     *           uy,vy,wy,uz,vz,wz                                      
#      data f0,f1,f3/0.d0,1.d0,3.d0/                                     
#      data pi2/6.283185307179586d0/                                     

    import numpy as np

    du = np.zeros(12)
    u = np.zeros(12)
                                                            
#======================================                                 
#=====  strike-slip contribution  =====                                 
#======================================                                 
    if(pot1 != f0):
        
        du[ 0] = alp1 * q / r3 + alp2 * x2 * qr                                
        du[ 1] = alp1 * x / r3 * sd + alp2 * xy * qr                                
        du[ 2] = -alp1 * x / r3 * cd + alp2 * x * d * qr                               
        du[ 3] = x * qr * (-alp1 + alp2 * (f1 + a5))                             
        du[ 4] = alp1 * a3 / r3 * sd + alp2 * y * qr * a5                             
        du[ 5] = -alp1 * a3 / r3 * cd + alp2 * d * qr * a5                             
        du[ 6] = alp1 * (sd / r3 - y * qr) + alp2 * f3 * x2 / r5 * uy                     
        du[ 7] = f3 * x / r5 * (-alp1 * y * sd + alp2 * (y * uy + q))                    
        du[ 8] = f3 * x / r5 * (alp1 * y * cd + alp2 * d * uy)                        
        du[ 9] = alp1 * (cd / r3 + d * qr) + alp2 * f3 * x2 / r5 * uz                     
        du[10] = f3 * x / r5 * (alp1 * d * sd + alp2 * y * uz)                        
        du[11] = f3 * x / r5 * (-alp1 * d * cd + alp2 * (d * uz - q))                    

        u = u + pot1 / pi2 * du

#===================================                                    
#=====  dip-slip contribution  =====                                    
#===================================                                    

    if(pot2 != f0):

        du[ 0] = alp2 * x * p * qr                                  
        du[ 1] = alp1 * s / r3 + alp2 * y * p * qr                                       
        du[ 2] = -alp1 * t / r3 + alp2 * d * p * qr                                          
        du[ 3] = alp2 * p * qr * a5                              
        du[ 4] = -alp1 * f3 * x * s / r5 - alp2 * y * p * qrx                                   
        du[ 5] = alp1 * f3 * x * t / r5 - alp2 * d * p * qrx                                    
        du[ 6] = alp2 * f3 * x / r5 * vy                           
        du[ 7] = alp1 * (s2d / r3 - f3 * y * s / r5) + alp2 * (f3 * y / r5 * vy + p * qr)                   
        du[ 8] = -alp1 * (c2d / r3 - f3 * y * t / r5) + alp2 * f3 * d / r5 * vy                       
        du[ 9] = alp2 * f3 * x / r5 * vz                           
        du[10] = alp1 * (c2d / r3 + f3 * d * s / r5) + alp2 * f3 * y / r5 * vz                       
        du[11] = alp1 * (s2d / r3 - f3 * d * t / r5) + alp2 * (f3 * d / r5 * vz - p * qr)                   

        u = u + pot2 / pi2 * du

#========================================                               
#=====  tensile-fault contribution  =====                               
#========================================                               
    if (pot3 != f0):
          
        du[ 0] = alp1 * x / r3 - alp2 * x * q * qr                              
        du[ 1] = alp1 * t / r3 - alp2 * y * q * qr                                   
        du[ 2] = alp1 * s / r3 - alp2 * d * q * qr                                      
        du[ 3] = alp1 * a3 / r3 - alp2 * q * qr * a5                          
        du[ 4] = -alp1 * f3 * x * t / r5 + alp2 * y * q * qrx                               
        du[ 5] = -alp1 * f3 * x * s / r5 + alp2 * d * q * qrx                                
        du[ 6] = -alp1 * f3 * xy / r5 - alp2 * x * qr * wy                          
        du[ 7] = alp1 * (c2d / r3 - f3 * y * t / r5) - alp2 * (y * wy + q) * qr                     
        du[ 8] = alp1 * (s2d / r3 - f3 * y * s / r5) - alp2 * d * qr * wy                      
        du[ 9] = alp1 * f3 * x * d / r5 - alp2 * x * qr * wz                          
        du[10] = -alp1 * (s2d / r3 - f3 * d * t / r5) - alp2 * y * qr * wz                      
        du[11] = alp1 * (c2d / r3 + f3 * d * s / r5) - alp2 * (d * wz - q) * qr                     

        u = u + pot3 / pi2 * du

#=========================================                              
#=====  inflate source contribution  =====                              
#=========================================                              
    if (pot4 != f0):
        du[ 0] = -alp1 * x / r3              
        du[ 1] = -alp1 * y / r3                      
        du[ 2] = -alp1 * d / r3                         
        du[ 3] = -alp1 * a3 / r3                  
        du[ 4] = alp1 * f3 * xy / r5                    
        du[ 5] = alp1 * f3 * x * d / r5                    
        du[ 6] = du[5]                                
        du[ 7] = -alp1 * b3 / r3                          
        du[ 8] = alp1 * f3 * y * d / r5                   
        du[ 9] = -du[6]                                
        du[10] = -du[9]                            
        du[11] = alp1 * c3 / r3                          

        u = u + pot4 / pi2 * du

    return u 


###############################################################################
def ub0(x, y, d, z, pot1, pot2, pot3, pot4):                    
###############################################################################
    import numpy as np

    u = np.zeros(12)
    du = np.zeros(12)

#       common /c0/alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d  
#       common /c1/p,q,s,t,xy,x2,y2,d2,r,r2,r3,r5,qr,qrx,a3,a5,b3,c3,     
#      *           uy,vy,wy,uz,vz,wz                                      
#       data f0,f1,f2,f3,f4,f5,f8,f9                                      
#      *        /0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,8.d0,9.d0/                 
#       data pi2/6.283185307179586d0/                                     

    c = d + z                                                             
    rd = r + d                                                            
    d12 = f1 / (r * rd * rd)                                                  
    d32 = d12 * (f2 * r + d) / r2                                               
    d33 = d12 * (f3 * r + d) / (r2 * rd)                                          
    d53 = d12 * (f8 * r2 + f9 * r * d + f3 * d2) / (r2 * r2 * rd)                           
    d54 = d12 * (f5 * r2 + f4 * r * d + d2) / r3 * d12

    fi1 = y * (d12 - x2 * d33)                                               
    fi2 = x * (d12 - y2 * d33)                                               
    fi3 = x / r3 - fi2                                                     
    fi4 = -xy * d32                                                       
    fi5 = f1 / (r * rd) - x2 * d32                                             
    fj1 = -f3 * xy * (d33 - x2 * d54)                                           
    fj2 = f1 / r3 - f3 * d12 + f3 * x2 * y2 * d54                                    
    fj3 = a3 / r3 - fj2                                                    
    fj4 = -f3 * xy / r5 - fj1                                                 
    fk1 = -y * (d32 - x2 * d53)                                               
    fk2 = -x * (d32 - y2 * d53)                                               
    fk3 = -f3 * x * d / r5 - fk2

    u = f0

#======================================                                 
#=====  strike-slip contribution  =====                                 
#======================================                                 
    if(pot1 != 0.):
        du[ 0] = -x2 * qr - alp3 * fi1 * sd                  
        du[ 1] = -xy * qr - alp3 * fi2 * sd                          
        du[ 2] = -c * x * qr - alp3 * fi4 * sd                             
        du[ 3] = -x * qr * (f1 + a5) - alp3 * fj1 * sd                 
        du[ 4] = -y * qr * a5 - alp3 * fj2 * sd                      
        du[ 5] = -c * qr * a5 - alp3 * fk1 * sd                       
        du[ 6] = -f3 * x2 / r5 * uy - alp3 * fj2 * sd                      
        du[ 7] = -f3 * xy / r5 * uy - x * qr - alp3 * fj4 * sd                     
        du[ 8] = -f3 * c * x / r5 * uy - alp3 * fk2 * sd                  
        du[ 9] = -f3 * x2 / r5 * uz + alp3 * fk1 * sd                          
        du[10] = -f3 * xy / r5 * uz + alp3 * fk2 * sd                      
        du[11] = f3 * x / r5 * (-c * uz + alp3 * y * sd)                        

        u = u + pot1 / pi2 * du                                        

#===================================                                    
#=====  dip-slip contribution  =====                                    
#===================================                                    
    if(pot2 != 0.):
        du[ 0] = -x * p * qr + alp3 * fi3 * sdcd              
        du[ 1] = -y * p * qr + alp3 * fi1 * sdcd                      
        du[ 2] = -c * p * qr + alp3 * fi5 * sdcd                         
        du[ 3] = -p * qr * a5 + alp3 * fj3 * sdcd                  
        du[ 4] = y * p * qrx + alp3 * fj1 * sdcd                       
        du[ 5] = c * p * qrx + alp3 * fk3 * sdcd                        
        du[ 6] = -f3 * x / r5 * vy + alp3 * fj1 * sdcd                   
        du[ 7] = -f3 * y / r5 * vy - p * qr + alp3 * fj2 * sdcd                  
        du[ 8] = -f3 * c / r5 * vy + alp3 * fk1 * sdcd               
        du[ 9] = -f3 * x / r5 * vz - alp3 * fk3 * sdcd                        
        du[10] = -f3 * y / r5 * vz - alp3 * fk1 * sdcd                    
        du[11] = -f3 * c / r5 * vz + alp3 * a3 / r3 * sdcd                     

        u = u + pot2 / pi2 * du                                        

#========================================                               
#=====  tensile-fault contribution  =====                               
#========================================                               
    if(pot3 != 0.):
        du[ 0] = x * q * qr - alp3 * fi3 * sdsd          
        du[ 1] = y * q * qr - alp3 * fi1 * sdsd                  
        du[ 2] = c * q * qr - alp3 * fi5 * sdsd                     
        du[ 3] = q * qr * a5 - alp3 * fj3 * sdsd              
        du[ 4] = -y * q * qrx - alp3 * fj1 * sdsd                   
        du[ 5] = -c * q * qrx - alp3 * fk3 * sdsd                    
        du[ 6] = x * qr * wy - alp3 * fj1 * sdsd                   
        du[ 7] = qr * (y * wy + q) - alp3 * fj2 * sdsd                  
        du[ 8] = c * qr * wy - alp3 * fk1 * sdsd               
        du[ 9] = x * qr * wz + alp3 * fk3 * sdsd                       
        du[10] = y * qr * wz + alp3 * fk1 * sdsd                   
        du[11] = c * qr * wz - alp3 * a3 / r3 * sdsd                    

        u = u + pot3 / pi2 * du

#=========================================                              
#=====  inflate source contribution  =====                              
#=========================================                              
    if (pot4 != 0.):
        du[ 0] = alp3 * x / r3       
        du[ 1] = alp3 * y / r3               
        du[ 2] = alp3 * d / r3                  
        du[ 3] = alp3 * a3 / r3           
        du[ 4] = -alp3 * f3 * xy / r5             
        du[ 5] = -alp3 * f3 * x * d / r5             
        du[ 6] = du[5]                         
        du[ 7] = alp3 * b3 / r3                   
        du[ 8] = -alp3 * f3 * y * d / r5            
        du[ 9] = -du[6]                         
        du[10] = -du[9]                     
        du[11] = -alp3 * c3 / r3                   
        
        u = u + pot4 / pi2 * du

    return(u)


###############################################################################
def uc0(x, y, d, z, pot1, pot2, pot3, pot4):                    
###############################################################################

    import numpy as np

    u = np.zeros(12)
    du = np.zeros(12)
                                                                     
#       common /c0/alp1,alp2,alp3,alp4,alp5,sd,cd,sdsd,cdcd,sdcd,s2d,c2d  
#       common /c1/p,q,s,t,xy,x2,y2,d2,r,r2,r3,r5,qr,qrx,a3,a5,b3,c3,um(6)
#       data f0,f1,f2,f3,f5,f7,f10,f15                                    
#      *        /0.d0,1.d0,2.d0,3.d0,5.d0,7.d0,10.d0,15.d0/               
#       data pi2/6.283185307179586d0/                                     

    c = d + z                                                             
    q2 = q * q                                                            
    r7 = r5 * r2                                                          
    a7 = f1 - f7 * x2 / r2                                                    
    b5 = f1 - f5 * y2 / r2                                                    
    b7 = f1 - f7 * y2 / r2                                                    
    c5 = f1 - f5 * d2 / r2                                                    
    c7 = f1 - f7 * d2 / r2                                                    
    d7 = f2 - f7 * q2 / r2                                                    
    qr5 = f5 * q / r2                                                       
    qr7 = f7 * q / r2                                                       
    dr5 = f5 * d / r2                                                       

#======================================                                 
#=====  strike-slip contribution  =====                                 
#======================================                                 
    if (pot1 != 0.):
        du[ 0] = -alp4 * a3 / r3 * cd + alp5 * c * qr * a5                           
        du[ 1] = f3 * x / r5 * (alp4 * y * cd + alp5 * c * (sd - y * qr5))                      
        du[ 2] = f3 * x / r5 * (-alp4 * y * sd + alp5 * c * (cd + d * qr5))                         
        du[ 3] = alp4 * f3 * x / r5 * (f2 + a5) * cd - alp5 * c * qrx * (f2 + a7)              
        du[ 4] = f3 / r5 * (alp4 * y * a5 * cd + alp5 * c * (a5 * sd - y * qr5 * a7))                 
        du[ 5] = f3 / r5 * (-alp4 * y * a5 * sd + alp5 * c * (a5 * cd + d * qr5 * a7))                  
        du[ 6] = du[4] / r5 * (alp4 * b5 * cd - alp5 * f5 * c / r2 * (f2 * y * sd + q * b7))                
        du[ 8] = f3 * x / r5 * (-alp4 * b5 * sd + alp5 * f5 * c / r2 * (d * b7 * sd - y * c7 * cd))          
        du[ 9] = f3 / r5 * (-alp4 * d * a5 * cd + alp5 * c * (a5 * cd + d * qr5 * a7))                  
        du[10] = f15 * x / r7 * (alp4 * y * d * cd + alp5 * c * (d * b7 * sd - y * c7 * cd))             
        du[11] = f15 * x / r7 * (-alp4 * y * d * sd + alp5 * c * (f2 * d * cd - q * c7))                   

        u = u + pot1 / pi2 * du
               
#===================================                                    
#=====  dip-slip contribution  =====                                    
#===================================                                    
    if (pot2 != 0.):
        du[ 0] = alp4 * f3 * x * t / r5 - alp5 * c * p * qrx                   
        du[ 1] = -alp4 / r3 * (c2d - f3 * y * t / r2) + alp5 * f3 * c / r5 * (s - y * p * qr5)         
        du[ 2] = -alp4 * a3 / r3 * sdcd + alp5 * f3 * c / r5 * (t + d * p * qr5)            
        du[ 3] = alp4 * f3 * t / r5 * a5 - alp5 * f5 * c * p * qr / r2 * a7      
        du[ 4] = f3 * x / r5 * (alp4 * (c2d - f5 * y * t / r2) - alp5 * f5 * c / r2 * (s - y * p * qr7))     
        du[ 5] = f3 * x / r5 * (alp4 * (f2 + a5) * sdcd - alp5 * f5 * c / r2 * (t + d * p * qr7))      
        du[ 6] = du[5]                                                           
        du[ 7] = f3 / r5 * (alp4 * (f2 * y * c2d + t * b5) + alp5 * c * (s2d - f10 * y * s / r2 - p * qr5 * b7))        
        du[ 8] = f3 / r5 * (alp4 * y * a5 * sdcd - alp5 * c * ((f3 + a5) * c2d + y * p * dr5 * qr7))     
        du[ 9] = f3 * x / r5 * (-alp4 * (s2d - t * dr5) - alp5 * f5 * c / r2 * (t + d * p * qr7))           
        du[10] = f3 / r5 * (-alp4 * (d * b5 * c2d + y * c5 * s2d) - alp5 * c * ((f3 + a5) * c2d + y * p * dr5 * qr7))    
        du[11] = f3 / r5 * (-alp4 * d * a5 * sdcd - alp5 * c * (s2d - f10 * d * t / r2 + p * qr5 * c7))       

        u = u + pot2 / pi2 * du

#========================================                               
#=====  tensile-fault contribution  =====                               
#========================================                               
    if(pot3 != 0.):
        du[ 0] = f3 * x / r5 * (-alp4 * s + alp5 * (c * q * qr5 - z))                     
        du[ 1] = alp4 / r3 * (s2d - f3 * y * s / r2) + alp5 * f3 / r5 * (c * (t - y + y * q * qr5) - y * z)
        du[ 2] = -alp4 / r3 * (f1 - a3 * sdsd) - alp5 * f3 / r5 * (c * (s - d + d * q * qr5) - d * z)
        du[ 3] = -alp4 * f3 * s / r5 * a5 + alp5 * (c * qr * qr5 * a7 - f3 * z / r5 * a5)          
        du[ 4] = f3 * x / r5 * (-alp4 * (s2d - f5 * y * s / r2) - alp5 * f5 / r2 * (c * (t - y + y * q * qr7) - y * z)) 
        du[ 5] = f3 * x / r5 * (alp4 * (f1 - (f2 + a5) * sdsd) + alp5 * f5 / r2 * (c * (s - d + d * q * qr7) - d * z)) 
        du[ 6] = du[5]                                                   
        du[ 7] = f3 / r5 * (-alp4 * (f2 * y * s2d + s * b5) - alp5 * (c * (f2 * sdsd + f10 * y * (t - y) / r2 - q * qr5 * b7) + z * b5))   
        du[ 8] = f3 / r5 * (alp4 * y * (f1 - a5 * sdsd) + alp5 * (c * (f3 + a5) * s2d - y * dr5 * (c * d7 + z)))             
        du[ 9] = f3 * x / r5 * (-alp4 * (c2d + s * dr5) + alp5 * (f5 * c / r2 * (s - d + d * q * qr7) - f1 - z * dr5))            
        du[10] = f3 / r5 * (alp4 * (d * b5 * s2d - y * c5 * c2d) + alp5 * (c * ((f3 + a5) * s2d - y * dr5 * d7) - y * (f1 + z * dr5)))        
        du[11] = f3 / r5 * (-alp4 * d * (f1 - a5 * sdsd) - alp5 * (c * (c2d + f10 * d * (s - d) / r2 - q * qr5 * c7) + z * (f1 + c5))) 

        u = u + pot3 / pi2 * du

#=========================================                              
#=====  inflate source contribution  =====                              
#=========================================                              
    if (pot4 != 0.):
        du[ 0] = alp4 * f3 * x * d / r5                
        du[ 1] = alp4 * f3 * y * d / r5                
        du[ 2] = alp4 * c3 / r3                    
        du[ 3] = alp4 * f3 * d / r5 * a5               
        du[ 4] = -alp4 * f15 * xy * d / r7              
        du[ 5] = -alp4 * f3 * x / r5 * c5               
        du[ 6] = du[5]                         
        du[ 7] = alp4 * f3 * d / r5 * b5               
        du[ 8] = -alp4 * f3 * y / r5 * c5               
        du[ 9] = du[6]                         
        du[10] = du[9]                         
        du[11] = alp4 * f3 * d / r5 * (f2 + c5)           

        u = u + pot4 / pi2 * du
        
    return(u)


###############################################################################
def dc3d0(alpha, \
###############################################################################
          x, y, z, \
          depth, \
          dip, \
          pot1, pot2, pot3, pot4, \
          ):

    """
        displacement and strain at depth                         
        due to buried point source in a semiinfinite medium      
                             coded by  y.okada ... sep.1991      
                             revised   y.okada ... nov.1991      

     input                                                            
       alpha : medium constant  (lambda+myu)/(lambda+2*myu)           
       x,y,z : coordinate of observing point                          
       depth : source depth                                           
       dip   : dip-angle (degree)                                     
       pot1-pot4 : strike-, dip-, tensile- and inflate-potency        
           potency=(  moment of double-couple  )/myu     for pot1,2   
           potency=(intensity of isotropic part)/lambda  for pot3     
           potency=(intensity of linear dipole )/myu     for pot4     
     output                                                           
       ux, uy, uz  : displacement ( unit=(unit of potency) /          
                   :                     (unit of x,y,z,depth)**2  )  
       uxx,uyx,uzx : x-derivative ( unit= unit of potency) /          
       uxy,uyy,uzy : y-derivative        (unit of x,y,z,depth)**3  )  
       uxz,uyz,uzz : z-derivative                                     
       iret        : return code  ( =0....normal,   =1....singular )  
    """

    global p, q, s, t, xy, x2, y2, d2, r, r2, r3, r5, qr, qrx, a3, a5, b3, c3, uy, vy, wy, uz, vz, wz
  
    import numpy as np
                                                              
#      common /c1/dummy(8),r,dumm(15)                                    
#      dimension  u(12),dua(12),dub(12),duc(12)                          

    if(z < 0.): print(" positive z was given in sub-dc3d0")
    
    u = np.zeros(12) + f0
    dua = np.zeros(12) + f0
    dub = np.zeros(12) + f0
    duc = np.zeros(12) + f0

    aalpha = alpha                                                      
    ddip = dip                                                          
    
    dccon0(aalpha, ddip)                                          

#======================================                                 
#=====  real-source contribution  =====                                 
#======================================                                 
    xx = x                                                              
    yy = y                                                              
    zz = z                                                              
    dd = depth + z                                                        
    
    dccon1(xx, yy, dd)                                             
    
    if (r == f0):
#=======================================                                
#=====  in case of singular (r=0)  =====                                
#=======================================                                
        ux = f0                                                             
        uy = f0                                                             
        uz = f0                                                             
        uxx = f0                                                            
        uyx = f0                                                            
        uzx = f0                                                            
        uxy = f0                                                            
        uyy = f0                                                            
        uzy = f0                                                            
        uxz = f0                                                            
        uyz = f0                                                            
        uzz = f0                                                            
        iret = 1   
                                                                 
        return(ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz)                                                            

#=======================================                                
#=====  case not singular (r != 0)  ====                                
#=======================================                                

#======================================
#=====  real-source contribution  =====
#======================================

    pp1 = pot1                                                          
    pp2 = pot2                                                          
    pp3 = pot3                                                          
    pp4 = pot4                                                          
    
    dua = ua0(xx, yy, dd, pp1, pp2, pp3, pp4)                            

    u[:9] = u[:9] - dua[:9] 
    u[9:] = u[9:] + dua[9:] 

#=======================================                                
#=====  image-source contribution  =====                                
#=======================================                                

    dd = depth - z                                                        
    dccon1(xx, yy, dd)                                             
    dua = ua0(xx, yy, dd, pp1, pp2, pp3, pp4)                            
    dub = ub0(xx, yy, dd, zz, pp1, pp2, pp3, pp4)                         
    duc = uc0(xx, yy, dd, zz, pp1, pp2, pp3, pp4)                         


    du = dua + dub + zz * duc                                      
    du[9:] = du[9:] + duc[0:3]
    u = u + du

    return u

