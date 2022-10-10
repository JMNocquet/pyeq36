class Dislocation:

    
    """
    Dislocation class
    Be careful in specifying coor_type 'xyz' means local cartesian (flat earth) , 'geo' means geographical coordinates
    x: Easting or longitude, y: Northing or latitude and depth is positive downward 
    """
###############################################################################
    def __init__ (self,index=None, \
                  x=None, y=None, depth=None, \
                  strike=None, dip=None, length=None, width=None, \
                  area=None, rake=None, max_slip=None, \
                  coor_type='xyz'):
###############################################################################
        
        self.index=index
        self.x=x
        self.y=y
        self.depth=depth
        self.strike=strike
        self.dip=dip
        self.length=length
        self.width=width
        self.area=area
        self.rake=rake
        self.max_slip=max_slip
        

###############################################################################
    def subfaults(self,n_length=None, n_width=None,coor_type='xyz',verbose=True):
###############################################################################
        """
        Divides a fault into subfaults and returns a list of Dislocation objects as defined in Dislocation 
        """
        
        import numpy as np
        
        if coor_type=='geo':
            deg2km=111.1
        else:
            deg2km=1.
        
        if verbose:
            np.set_printoptions(precision=2,linewidth=50   )
            print("-- Faults will be divided into ",n_length, "segments along strike and ",n_width, "segments along dip (",n_length*n_width," subfaults)")
        
        delta_length = self.length / float(n_length) / deg2km  
        delta_width =  self.width  / float(n_width)  / deg2km  
        
        print("-- delta_length & delta_width for calculation : %.3lf & %.3lf " % (delta_length,delta_width) )
        print("-- delta_length & delta_width for dislocation : %.3lf & %.3lf " % (delta_length*deg2km,delta_width*deg2km) )

        index=0
        lsubfaults=[]
        
        alpha = np.radians(-self.strike)
        
        Rdip = np.array([[np.cos(np.radians(self.dip)),0.,-np.sin(np.radians(self.dip))],\
                    [0.,1.,0.],\
                    [np.sin(np.radians(self.dip)),0.,np.cos(np.radians(self.dip))]])
        
        R2strike_hor = np.array([[np.cos(alpha),-np.sin(alpha),0.],\
                            [np.sin(alpha),np.cos(alpha),0.],[0.,0.,1.],])

        X0 = np.array([self.x,self.y,self.depth])
        
        for i in range(n_width):
            for j in range(n_length):
                
                local_x=i*delta_width
                local_y=j*delta_length
                
                LX = np.array([local_x,local_y,0.])

                X = np.dot(R2strike_hor, np.dot(Rdip,LX))

                X=X+X0
                
                
                depth = X0[2]+delta_width*deg2km*i*np.sin(np.radians(self.dip))
                
                length = delta_length*deg2km
                width = delta_width*deg2km
                area = length * width 
                
                if verbose:
                    print("-- subfault index : %04d depth: %.3lf length: %.3lf width %.3lf area: %.1lf " % (index,depth,length,width,area))
                
                subfault=Dislocation(index=index, x=X[0], y=X[1], \
                                     depth = depth, \
                                     strike=self.strike, dip=self.dip, \
                                     length = length, width = width , area = area )
                
                
                if verbose:
                    subfault.print_info()
                
                lsubfaults.append(subfault)
                index=index+1
                
        return(lsubfaults)
    
###############################################################################
    def corners(self,coor_type='xyz'):
###############################################################################
        """
        returns the 4 corners position vectors of a rectangular dislocation
        """

        from numpy import array
        from math import radians,cos,sin, pi
        
        if coor_type=='geo':deg2km=111.1
        else:deg2km=1.
        
        
        alpha=radians(-self.strike)+pi/2.
        
        # TOP 
        
        X1=array([self.x,self.y,self.depth])
        x2=self.x+cos(alpha)*self.length/deg2km
        y2=self.y+sin(alpha)*self.length/deg2km
        X2=array([x2,y2,self.depth])
        # BOTTOM
        delta_x=cos(radians(self.dip))*cos(alpha-pi/2.)*self.width/deg2km
        delta_y=cos(radians(self.dip))*sin(alpha-pi/2.)*self.width/deg2km
        delta_z=sin(radians(self.dip))*self.width
        DELTA=array([delta_x,delta_y,delta_z])
        
        X3=X2+DELTA
        X4=X1+DELTA
        return(X1,X2,X3,X4)
    
###############################################################################
    def centroid(self,coor_type='xyz'):
###############################################################################
        """
        Return the centroid of a rectangular dislocation
        """
        
        X1,X2,X3,X4=self.corners(coor_type=coor_type)
        return( (X1+X2+X3+X4)/4 )
    
###############################################################################
    def write_gmt(self,name,coor_type='xyz'):
###############################################################################
        """
        append rectangular fault information to a gmt psxy compatible file
        """
        
        import copy
        dis=copy.deepcopy(self)
        
        if coor_type=='geo':deg2km=111.1
        else:deg2km=1.

        (X1,X2,X3,X4)=self.corners()
        fs=open(name,'a+')
        fs.write('>\n')
        fs.write("%10.5lf  %10.5lf  \n"% (X1[0],X1[1]))
        fs.write("%10.5lf  %10.5lf  \n"%(X2[0],X2[1]))
        fs.write("%10.5lf  %10.5lf  \n"%(X3[0],X3[1]))
        fs.write("%10.5lf  %10.5lf  \n"%(X4[0],X4[1]))        
        fs.close()

###############################################################################
    def print_info(self):
###############################################################################
        """
        Print faults informations 
        """
        print(("index=%6d"%self.index))
        print(("(x,y,depth)=(%10.4lf,%10.4lf,%5.3lf )"%(self.x,self.y,self.depth)))
        print(("(strike,dip)=(%10.2lf,%10.2lf)" %(self.strike,self.dip)))
        print(("(length,width,area=(%5.2lf,%5.2lf,%5.2lf))"%(self.length,self.width,self.area)))
  
###############################################################################    
    def disp_slip_rake_no_edcmp(self,slip,rake,array_gps):
###############################################################################    
        
        # import
        import numpy as np
        from pyeq.green.okada_rde.okada import okada

        # elastic parameters        
        llambda = 3.2074E+10  
        mu      = 3.9701E+10

        # scaling for source point
        if self.area!=None and self.area > 0:
            true_area=self.area
            dislocation_area=self.width*self.length
            scaling=true_area / dislocation_area
            new_slip=slip*scaling
        else:
            new_slip=slip    

        # fault and slip parameters
        
        dislocations = np.zeros( 1 ) + new_slip
        xs = np.zeros( 1 ) + self.y
        ys = np.zeros( 1 ) + self.x
        zs = np.zeros( 1 ) + self.depth
        
        lengths = np.zeros( 1 ) + self.length
        widths  = np.zeros( 1 ) + self.width
        
        strikes = np.zeros( 1 ) + self.strike
        dips = np.zeros( 1 ) + self.dip
        
        if not rake:rake=self.rake

        rakes = np.zeros( 1 ) + rake
        
        # obs array parameters
        
        xrec = array_gps[:,1] 
        yrec = array_gps[:,0] 
        zrec0 = 0.

        ( disp, tilt , strain ) = okada(llambda, mu, dislocations, xs, ys, zs, lengths, widths, strikes, dips, rakes, xrec, yrec, zrec0)

        # reorder the output from edcmp
        # RR has the order x_east, y_north (in km) East, North, Up (positive upward) in millimeters
        # R  has the order x_north,y_east (in meters) D_North, D_East, D_Up (positive downward) in meters
        
        R = np.zeros( (xrec.shape[0], 5) )
        R[:,0] = xrec
        R[:,1] = yrec
        R[:,2:] = disp

        RR=np.zeros(R.shape)
        RR[:,0] =  R[:,1]*1.E-3
        RR[:,1] =  R[:,0]*1.E-3
        RR[:,2] =  R[:,3]
        RR[:,3] =  R[:,2]
        RR[:,4] = -R[:,4]
        return(RR)
  
###############################################################################
    def disp_slip_rake(self,slip,rake,array_gps):
###############################################################################
        """Calculates the 3D displacement along an array of points according to the
           displacement provided in slip and rake components.
           Array of points is a numpy array having x y
           Slip is rescaled by a factor equal to true_area/(length*width) if area is not None
        """
        
        import numpy as np
        import pyeq.lib.edcmp
        
        if not rake:rake=self.rake
        
        array_faults = np.array( [[ slip, self.x , self.y , self.depth,self.length,self.width,self.strike, self.dip,rake ]] )
        
        RR =  pyeq.lib.edcmp.disp_tilt_strain_stress_from_edcmp(array_faults,array_gps)

        return( RR )
        
        
        