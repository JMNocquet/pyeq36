"""
sequential solver with nnls
"""

def seq_nnls(ATA, ATB, nfaults, ntimestep , nnls , verbose):
    
    # m is the number of time
    
    # import 
    import numpy as np
    import pyeq.optimization.nnls.nnlsm
    from time import time

    my_init = np.zeros( nfaults )

    # loop on time step 
    for i in np.arange( ntimestep-1 )+1:
        print("-- step %04d / %04d" % (i,ntimestep) )
        
        sub_ATA = ATA[0:i*nfaults,0:i*nfaults]
        sub_ATB = ATB[0:i*nfaults]
        
        print("-- solving the linear system with algo")

        t0 = time()
        SLIP,_ = pyeq.optimization.nnls.nnlsm.nnlsm_blockpivot(sub_ATA, sub_ATB.reshape(-1, 1), is_input_prod=True, init=my_init.reshape(-1, 1))
        time_inversion=time()-t0
        print("-- time inversion blockpivot " , time_inversion )
        
        
#        t0 = time()
#        SLIP,_ = pyeq.lib.nnls.nnlsm.nnlsm_activeset(sub_ATA, sub_ATB.reshape(-1,1), is_input_prod=True, init=my_init.reshape(-1,1))
#        time_inversion=time()-t0
#        print("-- time inversion activeset " , time_inversion )
        
#        SLIP,time_inversion = pyeq.lib.make_inversion.pyeq_nnls(sub_ATA, sub_ATB, nnls, verbose=verbose) 
        

        my_init = np.append(SLIP,SLIP[-nfaults:]) 

        
        print("--saving slip for step %04d " % i)
        np.save('slip_'+str(i),SLIP)
