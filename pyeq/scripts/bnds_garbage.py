

# # CONSTRAINT ON THE CUMULATED SLIP AT A GIVEN TIME STEP
# 
# if args.c != None:
#     print '-- Dealing with absolute value constraints on slip'
#     for str_constraint in args.c:
#         index_subfault,cslip,time_step,slip_constraint=\
#         int(str_constraint.split('/')[0]),float(str_constraint.split('/')[1]),int(str_constraint.split('/')[2]),float(str_constraint.split('/')[3])
#         print("  -- Constraint on slip: subfault #%d at time step #%d : (%.1lf/%lf)" % (index_subfault,time_step,cslip,slip_constraint))
#         
#         G_CONSTRAINT=np.zeros((1,nstep*nfaults))
#         
#         for i in np.arange(np_dates.shape[0]-1):
#             if i == time_step:
#                 break
#             else:
#                 G_CONSTRAINT[0,i*SGEOMETRY.shape[0]+index_subfault]=1
# 
#         print np.sum(G_CONSTRAINT)
#         
#         d_CONSTRAINT=cslip
#         
#         G = np.append(G, G_CONSTRAINT, axis=0)
#         d = np.append(d, d_CONSTRAINT)
#         diag_Cd=np.append(diag_Cd,slip_constraint**2)
# 
# 
# # CONSTRAINT ON THE CUMULATED SLIP AT A GIVEN TIME STEP OVER A CIRCULAR AREA
# 
# if args.x != None:
#     print '-- Dealing with absolute value constraints on slip over a circular area'
#     for str_constraint in args.x:
#         
#         option=str_constraint.split('/')
#         
#         index_subfault,radius,cslip,time_step,constraint=map(float,option)
#         index_subfault=int(index_subfault)
# 
#         # find subfaults within area using Dm
#         lindex=np.argwhere(Dm[index_subfault]<radius).flatten()
#         
#         # loop on selected subfaults
#         
#         for index_subfault in lindex:
#         
#             print("  -- Constraint on slip: subfault #%d at time step #%d : (%.1lf/%lf)" % (index_subfault,time_step,cslip,constraint))
#         
#             G_CONSTRAINT=np.zeros((1,G.shape[1]))
#             
#             for i in np.arange(np_dates.shape[0]-1):
#                 G_CONSTRAINT[0,i*SGEOMETRY.shape[0]+index_subfault]=1
#                 if i == time_step:
#                     break
# 
#             d_CONSTRAINT=cslip
#             
#             G = np.append(G, G_CONSTRAINT, axis=0)
#             d = np.append(d, cslip)
#             diag_Cd=np.append(diag_Cd,constraint**2)
# 
# 
# 
# if args.e != None:
#     print '-- Dealing with equality constraints on slip evolution for two subfaults'
#     for str_constraint in args.e:
#         index_subfault1,index_subfault2,sigma_constraint=int(str_constraint.split('/')[0]),int(str_constraint.split('/')[1],float(str_constraint.split('/')[2]))
#         print("  -- Equality constraint for slip evolution at subfaults #%d = #%d (%lf)" % (index_subfault1,index_subfault2,sigma_constraint))
#         
#         
#         for i in np.arange(np_dates.shape[0]-1):
#             G_CONSTRAINT=np.zeros((1,G.shape[1]))
#             G_CONSTRAINT[0,i*SGEOMETRY.shape[0]+index_subfault1]=1.
#             G_CONSTRAINT[0,i*SGEOMETRY.shape[0]+index_subfault2]=-1.
#             d_CONSTRAINT=0.0
#         
#             G = np.append(G, G_CONSTRAINT, axis=0)
#             d = np.append(d, d_CONSTRAINT)
#             diag_Cd=np.append(diag_Cd,sigma_constraint**2)
# 
# 
# print '-- Dealing with the non-negativity constraints'
# 
# slip_no_constraint=100000.0 # 100 m
# 
# # SLIP_BOUNDS is a two column matrix having lower_bound and upper_bound in column 0 & 1
# 
# SLIP_BOUNDS=np.zeros((G.shape[1],2))
# SLIP_BOUNDS[:,1]=SLIP_BOUNDS[:,1]*0.0+slip_no_constraint
# print "    -- Slip bounds filled"
