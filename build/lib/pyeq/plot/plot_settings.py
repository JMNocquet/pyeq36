class plot_settings:
    def __init__ (self):
        # ncpu
        self.ncpu = 3
        # plot selection
        self.time_series = True
        self.cum_slip = False
        self.ccum_slip = True
        self.rate = False
        self.crate = True
        self.disp_obs = True
        self.disp_model = False
        self.disp_residual = False
        self.stf = True
        # plot options for time series
        self.ts_shift = False
        self.ts_min_yaxis = 10
        self.ts_ofmt = 'png'
        self.map = True
        # plot options for models
        self.shp_point = None
        self.shp_line = None
        self.shp_poly = None
        self.colormap = 'jet'
        self.ncolor = 20
        self.slip_bounds = None
        self.rate_bounds = None
        self.disp_scale = 1
        self.ncontour = 2
        self.interpolate = 'cubic'

    #    @classmethod
    ###################################################################
    def load(cls, file_name=None, verbose=False):
    ###################################################################
        """
        load an existing configuration for plot.

        param file_name: model file as text file
        param verbose: verbose mode

        """

        import pyeq.message.message as MESSAGE
        import pyeq.message.verbose_message as VERBOSE
        import pyeq.message.error as ERROR
        import pyeq.message.warning as WARNING
        import pyeq.message.debug_message as DEBUG

        try:
            cf = open(file_name,'r')
        except:
            ERROR("Could not open %s" % file_name)

        for line in cf:
            ###############################################################################
            # not to be read
            ###############################################################################

            # comment
            if line[0] == '#': continue
            # blank line
            if len(line.strip()) == 0: continue
            # incomplete line
            lline = line.split('=')
            if len(lline) < 2:
                ERROR(("reading file: %s line |%s|" % (file_name, line)), exit=True)

            key   = lline[0].strip()

            # get the value
            # get the position of the first occurrence of =
            try:
                idx = line.index('=')
            except:
                ERROR(("= (equal) character missing for line: %s " % line), exit=True)

            value = line[idx + 1:].strip()
            DEBUG(("reading key in conf : %s" % key))

            if key not in cls.__dict__:
                ERROR("plot option %s not recognized." % key)
            else:
                # case for shp_line and shp_poly which can be repeated
                if key in ['shp_line','shp_poly','shp_point']:
                    DEBUG("shp case %s,%s" % (key,value) )
                    DEBUG("cls.__dict__[key] %s" % cls.__dict__[key])
                    DEBUG(cls.__dict__[key] is None )

                    if cls.__dict__[key] in ['None','NONE'] or cls.__dict__[key] is None:
                        DEBUG("setting shp as []")
                        cls.__dict__[key] = []
                    if value not in ['None','NONE']:
                        cls.__dict__[key].append( value )
                    else:
                        pass
                else:
                    cls.__dict__[key] = value
#        except:
#            ERROR(("reading %s" % file_name), exit=True)

    def print(cls):
        """
        print current plot_settings object content
        """
        print("# PYEQ PLOT SETTING DEFAULTS")
        for k,v in cls.__dict__.items():
            print("%-15s = %s" %(k,v))