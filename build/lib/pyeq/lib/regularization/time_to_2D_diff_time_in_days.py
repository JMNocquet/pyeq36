def time_to_2D_diff_time_in_days( delta_time_in_days ):
    """
    Given a vector of dates in days, compute a matrix of time distance between each date pair 
    """

    import numpy as np
    
    
    delta_day_2D = np.zeros( ( delta_time_in_days.shape[0] , delta_time_in_days.shape[0]  ) )


    for i in np.arange( delta_day_2D.shape[0] ):
        for j in np.arange( i , delta_day_2D.shape[0] ):
            delta_day_2D[i,j] = delta_time_in_days[j] - delta_time_in_days[i]
            delta_day_2D[j,i] = delta_day_2D[i,j] 

    return delta_day_2D
