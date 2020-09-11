def display_array( array , comment = ''):
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    print('DISPLAY')
    plt.imshow(array, cmap="jet", aspect=1)
    plt.colorbar()
    title = ( "%s      %d x %d     max = %.2lf "  %  ( comment , array.shape[0] , array.shape[1] , np.max(array)      ) ) 
    plt.title( title )
    plt.show()