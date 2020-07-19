import numpy as np

# Helpers
def Expand_Contact_Matrices(contactMatrix, p_home, p_reduced, p_full):
    ''' 
    This function expands the 3x3 contact matrices into 5x5 ones by separating
    adult contacts based on home, reduced, or full contact probabilities.
    '''
    
    temp_rowexpand = contactMatrix[(0,1,1,1,2),:]
    temp_expand = temp_rowexpand[:,(0,1,1,1,2)]
    
    temp_expand[:,(1,2,3)] = temp_expand[:,(1,2,3)] * np.array((p_home, p_reduced, p_full)).reshape((1,3))
    
    return(temp_expand)

def Expand_10x10(MatA, MatB):
    '''
    This function interlaces a 5x5 high- and 5x5 low- contact matrix into a single 10x10
    '''
    ret = np.zeros((5,10))
    
    # fill columns
    ret[:,0::2] = MatA
    ret[:,1::2] = MatB
    
    # fill rows
    ret = ret[np.repeat(np.arange(5), 2),:]
    
    return(ret)