class pyeq_model:     
    def __init__ (self):
        pass

    @classmethod
    ###################################################################
    def load(cls, file_name=None, verbose=False):
    ###################################################################
        """
        load an existing pyeq model.

        param file_name: model file as pickle
        param verbose: verbose mode

        """

        import pickle
        import os
        from colors import red
        import sys

        import pyeq.message.message as MESSAGE
        import pyeq.message.verbose_message as VERBOSE
        import pyeq.message.error as ERROR
        import pyeq.message.warning as WARNING
        import pyeq.message.debug_message as DEBUG

        try:
            MESSAGE("Loading %s (%.2f Gb) " % ( file_name , os.path.getsize( file_name ) /1024 / 1024 / 1024 ) )
            with open( file_name, "rb") as f:
                model = pickle.load( f )
            f.close()
            VERBOSE("model object loaded.")
        except:
            ERROR(("Could not load: %s " % ( file_name ) ),exit=True)

        return model

    ###################################################################
    def print_info(cls):
    ###################################################################
        """
        print info on a model
        """

        import numpy as np

        import pyeq.message.message as MESSAGE
        import pyeq.message.verbose_message as VERBOSE
        import pyeq.message.error as ERROR
        import pyeq.message.warning as WARNING
        import pyeq.message.debug_message as DEBUG

        for attr in sorted( cls.__dict__):
            if isinstance( cls.__dict__[attr], float ):
                MESSAGE("%s : %.2E" % (attr,cls.__dict__[attr]))
            if isinstance( cls.__dict__[attr], str ):
                MESSAGE("%s : %s" % (attr,cls.__dict__[attr]))
            if isinstance( cls.__dict__[attr], np.ndarray ):
                MESSAGE("%s : %s"  % (attr,  cls.__dict__[attr].shape,))

    ###################################################################
    def write_pickle(cls, file_name=None):
    ###################################################################
        """
        write a model object as pickle
        """

        import pickle
        import os

        import pyeq.message.message as MESSAGE
        import pyeq.message.verbose_message as VERBOSE
        import pyeq.message.error as ERROR
        import pyeq.message.warning as WARNING
        import pyeq.message.debug_message as DEBUG

        VERBOSE("writing %s" % file_name  )

        ofile = open( file_name , 'wb')
        pickle.dump( cls , ofile , pickle.DEFAULT_PROTOCOL)
        ofile.close()

        VERBOSE("%s (%.2f Gb) written" % (file_name, os.path.getsize(file_name) / 1024 / 1024 / 1024))

