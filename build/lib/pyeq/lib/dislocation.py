class Dislocation:

    
    """Dislocation class: creates an elastic dislocation 
    """

###############################################################################    
    def __init__ (self,index=None, x=None, y=None, depth=None, \
###############################################################################    
                  strike=None, dip=None, length=None, width=None,\
                   area=None, rake=None, max_slip=None,\
                   centroid_x=None,centroid_y=None,centroid_z=None):
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
        self.centroid_x=centroid_x
        self.centroid_y=centroid_y
        self.centroid_z=centroid_z
        
    
###############################################################################    
    def centroid_to_upper_left_corner(self):
###############################################################################    
        """
        From the centroid coordinates length, width, strike and dip, returns the coordinates of the upper left corner
        """
        
        from numpy import array
        from math import radians,cos,sin, pi
        alpha=radians(-self.strike)+pi/2.
        
        # CENTROID
        
        CENTROID=array([self.centroid_x,self.centroid_y,self.centroid_z])
        
        # ALONG_STRIKE
        
        delta_x=cos(alpha)*self.length
        delta_y=sin(alpha)*self.length
        delta_z=0.0
        DELTA_ALONG_STRIKE=array([delta_x,delta_y,delta_z])
        
        # ALONG DIP
        delta_x=cos(radians(self.dip))*cos(alpha-pi/2.)*self.width
        delta_y=cos(radians(self.dip))*sin(alpha-pi/2.)*self.width
        delta_z=sin(radians(self.dip))*self.width
        DELTA_ALONG_DIP=array([delta_x,delta_y,delta_z])
        
        TOP_LEFT_CORNER=CENTROID-DELTA_ALONG_STRIKE/2.-DELTA_ALONG_DIP/2.
        
        self.x=TOP_LEFT_CORNER[0]
        self.y=TOP_LEFT_CORNER[1]
        self.depth=TOP_LEFT_CORNER[2]
        
        return(TOP_LEFT_CORNER)

        
        

###############################################################################    
    def subfaults(self,n_length=None, n_width=None):
###############################################################################    
        """Divides a rectangular fault into subfaults and returns a list of Dislocation objects as defined in Dislocation 
        """
        from numpy import array,dot,set_printoptions
        set_printoptions(precision=2,linewidth=50   )
        from math import cos,sin,radians
print("-- Faults will be divided into ",n_length, "segments along strike and ",n_width, "segments along dip (",n_length*n_width," subfaults)")
        
        delta_length=self.length/float(n_length)
        delta_width=self.width/float(n_width)
        index=0
        lsubfaults=[]
        alpha=radians(-self.strike)
        Rdip=array([[cos(radians(self.dip)),0.,-sin(radians(self.dip))],[0.,1.,0.],[sin(radians(self.dip)),0.,cos(radians(self.dip))]])
        #print 'Rdip ',Rdip
        R2strike_hor=array([[cos(alpha),-sin(alpha),0.],[sin(alpha),cos(alpha),0.],[0.,0.,1.],])
        #print 'R2strike_hor ',R2strike_hor
        R=dot(Rdip,R2strike_hor)
        #print 'R ',R

        X0=array([self.x,self.y,self.depth])
        #print 'X0 ',X0
        for i in range(n_width):
            #print 'i ',i
            for j in range(n_length):
                #print 'j ',j
                # faults coordinates
                local_x=i*delta_width
                local_y=j*delta_length
                LX=array([local_x,local_y,0.])
                #print 'LX ',LX
                #print 'RDip X ',dot(Rdip,LX)
                X=dot(R2strike_hor,dot(Rdip,LX))
                #print 'dX ',X
                X=X+X0
                #print 'X ',X
                
                #print "local x,y, X",local_x,local_y,X+X0
                
                subfault=Dislocation(index=index, x=X[0], y=X[1], depth=X[2], strike=self.strike, dip=self.dip, length=delta_length, width=delta_width)
                
                lsubfaults.append(subfault)
                
                
                
                index=index+1
                
                
                
        return(lsubfaults)
    
###############################################################################    
    def corners(self):
###############################################################################    
        """returns the 4 corners position vectors"""
        from numpy import array
        from math import radians,cos,sin, pi
        alpha=radians(-self.strike)+pi/2.
        
        # TOP
        
        X1=array([self.x,self.y,self.depth])
        x2=self.x+cos(alpha)*self.length
        y2=self.y+sin(alpha)*self.length
        X2=array([x2,y2,self.depth])
        # BOTTOM
        delta_x=cos(radians(self.dip))*cos(alpha-pi/2.)*self.width
        delta_y=cos(radians(self.dip))*sin(alpha-pi/2.)*self.width
        delta_z=sin(radians(self.dip))*self.width
        DELTA=array([delta_x,delta_y,delta_z])
        
        X3=X2+DELTA
        X4=X1+DELTA
        return(X1,X2,X3,X4)
    
###############################################################################    
    def centroid(self):
###############################################################################    
        """
        Returns the dislocation centroid
        """
        X1,X2,X3,X4=self.corners()
        X=(X1+X2+X3+X4)/4.
        
        self.centroid_x=X[0]
        self.centroid_y=X[1]
        self.centroid_z=X[2]
        
        return(X)
        
    
###############################################################################    
    def write_gmt(self,name):
###############################################################################    
        """write fault information to a gmt psxy compatible file"""
        (X1,X2,X3,X4)=self.corners()
        fs=open(name,'a+')
        fs.write('>\n')
        xe=self.x
        fs.write("%10.5lf  %10.5lf  \n"% (X1[0],X1[1]))
        fs.write("%10.5lf  %10.5lf  \n"%(X2[0],X2[1]))
        fs.write("%10.5lf  %10.5lf  \n"%(X3[0],X3[1]))
        fs.write("%10.5lf  %10.5lf  \n"%(X4[0],X4[1]))        
        fs.close()

###############################################################################    
    def print_info(self):
###############################################################################    
        """Print faults informations 
        """
        print(("index=%6d"%self.index))
        print(("(x,y,depth)=(%10.4lf,%10.4lf,%5.3lf )"%(self.x,self.y,self.depth)))
        print(("(strike,dip)=(%10.2lf,%10.2lf)" %(self.strike,self.dip)))
        print(("(length,width=(%5.2lf,%5.2lf))"%(self.length,self.width)))

###############################################################################    
    def disp_slip_rake_no_edcmp(self,slip,rake,array_gps):
###############################################################################    
        
        # import
        import numpy as np
        
        
        from pyeq.green.okada_rde.okada import okada
        
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
        """Calculates the Green function for the dislocation, i.e. the 3D displacement along an array of points according to the
           displacement provided in slip and rake components.
           Array of points is a numpy array having x y
           Slip is rescaled by a factor equal to true_area/(length*width)
           uses the edcmp
           The return is
        """
        
        if not rake:rake=self.rake
        # first prepare edcmp input file
           
        fedcmp_input='input_edcmp.dat'
        fs=open(fedcmp_input,'w')
           
        intro_edcmp =''
        intro_edcmp += '#==============================================================================='+"\n"
        intro_edcmp += '# This is the input file of FORTRAN77 program "edcomp" for calculating'+"\n"
        intro_edcmp += '# earthquakes static deformations (3 displacement components, 6 strain/stress'+"\n"
        intro_edcmp += '# components and 2 vertical tilt components) based on the dislocation theory.'+"\n"
        intro_edcmp += '# The earth model used is either a homogeneous or multi-layered, isotropic and'+"\n"
        intro_edcmp += '# elastic half-space. The earthquake source is represented by an arbitrary number'+"\n"
        intro_edcmp += '# of rectangular dislocation planes.'+"\n"
        intro_edcmp += '#'+"\n"
        intro_edcmp += '# Note that the Cartesian coordinate system is used and the seismological'+"\n"
        intro_edcmp += '# convention is adopted, that is, x is northward, y is eastward, and z is'+"\n"
        intro_edcmp += '# downward.'+"\n"
        intro_edcmp += '#'+"\n"
        intro_edcmp += '# First implemented in Potsdam, Feb, 1999'+"\n"
        intro_edcmp += '# Last modified: Potsdam, Nov, 2001'+"\n"
        intro_edcmp += '#'+"\n"
        intro_edcmp += '# by'+"\n"
        intro_edcmp += '# Rongjiang Wang, Frank Roth, & Francisco Lorenzo'+"\n"
        intro_edcmp += '# GeoForschungsZetrum Potsdam, Telegrafenberg, 14473 Potsdam, Germany'+"\n"
        intro_edcmp += '#'+"\n"
        intro_edcmp += '# For questions and suggestions please send e-mails to wang@gfz-potsdam.de'+"\n"
        intro_edcmp += '#==============================================================================='+"\n"

        obs_edcmp =''
        obs_edcmp += '# OBSERVATION ARRAY'+"\n"
        obs_edcmp += '# ================='+"\n"
        obs_edcmp += '# 1. switch for irregular positions (0) or a 1D profile (1)'+"\n"
        obs_edcmp += '#    or a rectangular 2D observation array (2): ixyr'+"\n"
        obs_edcmp += '#'+"\n"
        obs_edcmp += '#    IF (1 for irregular observation positions) THEN'+"\n"
        obs_edcmp += '#    '+"\n"
        obs_edcmp += '# 2. number of positions: nr'+"\n"
        obs_edcmp += '# 3. coordinates of the observations: (xr(i),yr(i)),i=1,nr'+"\n"
        obs_edcmp += '#'+"\n"
        obs_edcmp += '#    ELSE IF (the switch = 1 for a 1D profile) THEN'+"\n"
        obs_edcmp += '#'+"\n"
        obs_edcmp += '# 2. number of position samples of the profile: nr'+"\n"
        obs_edcmp += '# 3. the start and end positions: (xr1,yr1), (xr2,yr2)'+"\n"
        obs_edcmp += '#'+"\n"
        obs_edcmp += '#    ELSE IF (2 for rectanglular 2D observation array) THEN'+"\n"
        obs_edcmp += '#'+"\n"
        obs_edcmp += '# 2. number of xr samples, start and end values [m]: nxr, xr1,xr2'+"\n"
        obs_edcmp += '# 3. number of yr samples, start and end values [m]: nyr, yr1,yr2'+"\n"
        obs_edcmp += '#'+"\n"
        obs_edcmp += '#    Note that the total number of observation positions (nr or nxr*nyr)'+"\n"
        obs_edcmp += '#    should be <= NRECMAX (see edcglobal.h)!'+"\n"
        obs_edcmp += '#==============================================================================='+"\n"
        obs_edcmp += '#'+"\n"

        output_edcmp=''
        output_edcmp += '#==============================================================================='+"\n"
        output_edcmp += '# OUTPUTS'+"\n"
        output_edcmp += '# ======='+"\n"
        output_edcmp += '# 1. output directory in char format: outdir'+"\n"
        output_edcmp += '# 2. select the desired outputs (1/0 = yes/no)'+"\n"
        output_edcmp += '# 3. the file names in char format for displacement vector, strain tensor,'+"\n"
        output_edcmp += '#    stress tensor, and vertical tilts:'+"\n"
        output_edcmp += '#    dispfile, strainfile, stressfile, tiltfile'+"\n"
        output_edcmp += '#'+"\n"
        output_edcmp += '#    Note that all file or directory names should not be longer than 80'+"\n"
        output_edcmp += '#    characters. Directories must be ended by / (unix) or \ (dos)!'+"\n"
        output_edcmp += '#==============================================================================='+"\n"
        output_edcmp += '  \'./\''+"\n"
        output_edcmp += '        1               0              0              0'+"\n"

        source_edcmp =''

        source_edcmp += '#==============================================================================='+"\n"
        source_edcmp += '# RECTANGLAR DISLOCATION SOURCES'+"\n"
        source_edcmp += '# =============================='+"\n"
        source_edcmp += '# 1. number of the source rectangles: ns (<= NSMAX in edcglobal.h)'+"\n"
        source_edcmp += '# 2. the 6 parameters for the 1. source rectangle:'+"\n"
        source_edcmp += '#    Slip [m],'+"\n"
        source_edcmp += '#    coordinates of the upper reference point for strike (xs, ys, zs) [m],'+"\n"
        source_edcmp += '#    length (strike direction) [m], and width (dip direction) [m],'+"\n"
        source_edcmp += '#    strike [deg], dip [deg], and rake [deg]'+"\n"
        source_edcmp += '# 3. ... for the 2. source ...'+"\n"
        source_edcmp += '# ...'+"\n"
        source_edcmp += '#                   N'+"\n"
        source_edcmp += '#                  /'+"\n"
        source_edcmp += '#                 /| strike'+"\n"
        source_edcmp += '#         Ref:-> @------------------------'+"\n"
        source_edcmp += '#                |\        p .            \ W'+"\n"
        source_edcmp += '#                :-\      i .              \ i'+"\n"
        source_edcmp += '#                |  \    l .                \ d'+"\n"
        source_edcmp += '#                :90 \  S .                  \ t'+"\n"
        source_edcmp += '#                |-dip\  .                    \ h'+"\n"
        source_edcmp += '#                :     \. | rake               \ '+"\n"
        source_edcmp += '#                Z      -------------------------'+"\n"
        source_edcmp += '#                              L e n g t h'+"\n"
        source_edcmp += '#'+"\n"
        source_edcmp += '#    Note that if one of the parameters length and width = 0, then a line source'+"\n"
        source_edcmp += '#    will be considered and the displocation parameter Slip has the unit m^2 if'+"\n"
        source_edcmp += '#    both length and width = 0, then a point source will be considered and the'+"\n"
        source_edcmp += '#    Slip has the unit m^3.'+"\n"
        source_edcmp += '#==============================================================================='+"\n"


        fault_edcmp = ''
        fault_edcmp += '#         coord. origin: (40.739N, 30.05E)'+"\n"
        fault_edcmp += '#-------------------------------------------------------------------------------'+"\n"
        fault_edcmp += '# no  Slip   xs        ys       zs        length    width   strike   dip  rake'+"\n"
        fault_edcmp += '#-------------------------------------------------------------------------------'+"\n"
        fault_edcmp += '#   1   1.00 0.0d+0  50.0d+03  0.0d+03   64.0d+03  15.0d+03  45.0   45.0  -90.0'+"\n"

        earth_model_edcmp=''
        earth_model_edcmp +=  '#==============================================================================='+"\n"
        earth_model_edcmp += '# If the earth model used is a layered half-space, then the numerical Green\'s'+"\n"
        earth_model_edcmp += '# function approach is applied. The Green\'s functions should have been prepared'+"\n"
        earth_model_edcmp += '# with the program "edgrn" before the program "edcmp" is started. In this case,'+"\n"
        earth_model_edcmp += '# the following input data give the addresses where the Green\'s functions have'+"\n"
        earth_model_edcmp += '# been stored and the grid side to be used for the automatic discretization'+"\n"
        earth_model_edcmp += '# of the finite rectangular sources.'+"\n"
        earth_model_edcmp += '#'+"\n"
        earth_model_edcmp += '# If the earth model used is a homogeneous half-space, then the analytical'+"\n"
        earth_model_edcmp += '# method of Okada (1992) is applied. In this case, the Green\'s functions are'+"\n"
        earth_model_edcmp += '# not needed, and the following input data give the shear modulus and the'+"\n"
        earth_model_edcmp += '# Poisson ratio of the model.'+"\n"
        earth_model_edcmp += '#==============================================================================='+"\n"
        earth_model_edcmp += '# CHOICE OF EARTH MODEL'+"\n"
        earth_model_edcmp += '# ====================='+"\n"
        earth_model_edcmp += '# 1. switch for layered (1) or homogeneous (0) model'+"\n"
        earth_model_edcmp += '#'+"\n"
        earth_model_edcmp += '#    IF (layered model) THEN'+"\n"
        earth_model_edcmp += '#'+"\n"
        earth_model_edcmp += '# 2. directory of the Green\'s functions and the three files for the'+"\n"
        earth_model_edcmp += '#    fundamental Green\'s functions: grndir, grnfiles(3)'+"\n"
        earth_model_edcmp += '#'+"\n"
        earth_model_edcmp += '#    Note that all file or directory names should not be longer than 80'+"\n"
        earth_model_edcmp += '#    characters. Directories must be ended by / (unix) or \ (dos)!'+"\n"
        earth_model_edcmp += '#'+"\n"
        earth_model_edcmp += '#    ELSE (homogeneous model) THEN'+"\n"
        earth_model_edcmp += '#'+"\n"
        earth_model_edcmp += '# 2. the observation depth, the two Lame constants parameters of the homogeneous'+"\n"
        earth_model_edcmp += '#    model: zrec [m], lambda [Pa], mu [Pa]'+"\n"
        earth_model_edcmp += '#==============================================================================='+"\n"
        earth_model_edcmp += '#  1'+"\n"
        earth_model_edcmp += '#  \'./grnfcts/\'  \'izmhs.ss\'  \'izmhs.ds\'  \'izmhs.cl\''+"\n"
        earth_model_edcmp += '  0'+"\n"
        earth_model_edcmp += '  0.00d+00  0.28758E+11  0.29353E+11'+"\n"
        earth_model_edcmp += '#================================end of input==================================='+"\n"
        
        fs.write(intro_edcmp)
        fs.write(obs_edcmp)
        fs.write("  0\n")
        # edcmp convention X = north, Y = East
        fs.write("  %d\n" % array_gps.shape[0])
        line=0
        for i in range(array_gps.shape[0]-1):
            #print 'i ',i
            #print 'array_gps ',array_gps
            fs.write("  (%5.4lfd+03,%5.4lfd+03)," %(array_gps[i,1],array_gps[i,0]))
            line+=1
            if (line == 3):
                line=0;fs.write("\n")
        i=array_gps.shape[0]-1        
        fs.write("  (%5.4lfd+03,%5.4lfd+03)" %(array_gps[i,1],array_gps[i,0]))                
        fs.write("\n")
        fs.write(output_edcmp)
        fs.write("  \'edcmp.disp\'    \'edcmp.strn\'   \'edcmp.strss\'  \'edcmp.tilt\'"+"\n")
        fs.write(source_edcmp)
        fs.write(" 1\n")
        fs.write(fault_edcmp)
        # scaling for source point
        if self.area!=None and self.area > 0:
            true_area=self.area
            dislocation_area=self.width*self.length
            scaling=true_area / dislocation_area
            new_slip=slip*scaling
        else:
            new_slip=slip    
        
        fs.write("  1   %5.2lf %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lf %5.1lf %5.1lf\n" %(new_slip,self.y,self.x,self.depth,self.length,self.width,self.strike,self.dip,rake))
        fs.write(earth_model_edcmp)
        fs.close()
        
        # now runs edcmp
        fs=open('linput_edcmp','w')
        fs.write("%s\n" % fedcmp_input)
        fs.close()
        #from subprocess import call
        #call("edcmp < linput_edcmp", shell=False) 
        #call(['edcmp < linput_edcmp'])
        import subprocess
        subprocess.getstatusoutput('edcmp < linput_edcmp')
        # reads results
        import numpy as np
        R=np.genfromtxt('edcmp.disp',dtype='float', comments='#')
        
        #### if edcmp.disp only has one record a 1-D array (vector) is returned to R instead of a 2-D array
        if (R.ndim ==1): R.resize(1,5)
        
        # reorder the output from edcmp
        # RR has the order x_east, y_north (in km) East, North, Up (positive upward) in millimeters
        # R  has the order x_north,y_east (in meters) D_North, D_East, D_Up (positive downward) in meters
        RR=np.zeros(R.shape)
        RR[:,0]=R[:,1]*1.E-3
        RR[:,1]=R[:,0]*1.E-3
        RR[:,2]=R[:,3]
        RR[:,3]=R[:,2]
        RR[:,4]=-R[:,4]
        #print RR
        return(RR)

        
        