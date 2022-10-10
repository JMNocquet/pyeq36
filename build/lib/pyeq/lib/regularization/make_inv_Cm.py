#def make_inv_Cm_step(model ):
#
#
#    if model.Cm_normalize == 'd0':
#
#        # gets dc0 = minimum characteristic length for fault discretization
#        dc0 = np.sqrt(np.min( model.SGEOMETRY.rdis_area))
#    else:
#        dc0 = 1.
#
#
#    # case
#    if model.Cm_type == 'exponential':
#        CM=np.exp(-Dm/model.dc) * model.sdotmax**2 * w
#    if model.Cm_type == 'm_exponential':
#        CM=np.exp(-Dm/model.dc) * model.sdotmax**2 * w * (dc0/model.dc)**2
#
#    inv_CM_step = pyacs.lib.glinalg.cov_to_invcov(CM)
        