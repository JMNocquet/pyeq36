"""
Package that allow to create and manipulate the Green tensor relating slip to displacement, strain and stress.
By default, GREEN is a TENSOR OF DIM 4:
    GREEN(i,j,k,l) is the prediction for dislocation i at site j component k for rake l
    k=0,1,2 = east, north, up
    l=0,1 : rake_00 & rake_90; optionally l=2 for tensile slip

Current implementation includes:

- Okada for rectangular dislocations in a homogeneous half-space
- edcmp for rectangular dislocations in a homogeneous half-space
- Meade for triangular dislocations in a homogeneous half-space
- Nikkhoo for triangular dislocations in a homogeneous half-space
- Nikkhoo for rectangular dislocations in a homogeneous half-space (from assembling two triangular dislocations)
"""
