def print_model_attributes( model ):
    
    for att in sorted( dir( model )):
        if att[0] != '_':
            print("model.%s" % att)

def print_diff_model_attributes( model_old , model_new ):
    """
    print model attributes differences
    """
    
    import numpy as np

    def isRecarray(a):
        return a.dtype.fields != None
    
    latt_old = []
    for att in model_old.__dict__:
        latt_old.append( att )
            
    latt_new = []
    for att in model_new.__dict__:
        latt_new.append( att )
    
    # print new attributes
    
    latt_common = []
    
    for att in latt_new:
        if att not in latt_old:
            print("-- attribute %25s     added" % att )
        else:
            latt_common.append( att )
    
    # print delete attribute

    for att in latt_old:
        if att not in latt_new:
            print("-- attribute %25s     deleted" % att )

    # print changed attributes

    for att in latt_common:
        # numpy array case (not recarray)
        if not hasattr(model_old,att):
            print("-- attribute %25s     missing" % att )
            next()

        if isinstance( model_old.__dict__[ att ], np.ndarray):
            if isRecarray( model_old.__dict__[ att ] ):
                continue
                
            # check their shape
            if model_old.__dict__[ att ].shape == model_new.__dict__[ att ].shape:
                
                # check they are string
                if model_old.__dict__[ att ].dtype.kind == 'U':
                    # test for string array
                    
                    if not np.array_equal(  model_old.__dict__[ att ] , model_new.__dict__[ att ]  ):
                         print("-- attribute numpy array %14s     updated" % att )
                else:
                    # test for numeric array
                    if not np.allclose( model_old.__dict__[ att ] , model_new.__dict__[ att ] , equal_nan=True ):
                        print("-- attribute numpy array %14s     updated" % att )
            else:
                # not the same size so att was changed
                print("-- attribute numpy array %14s     updated" % att )
        else:
            # not hte same type so att was updated
            if model_old.__dict__[ att ] != model_new.__dict__[ att ]:
                print("-- attribute %25s     updated" % att )
    
     
    
    