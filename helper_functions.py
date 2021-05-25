def UPGMA(distances):
    """Unweighted pair group method with arithmetic mean (UPGMA) agglomerative clustering.
    
    Parameters
    ----------
    distances: np.ndarray
        A two dimensional, square, symmetric matrix containing distances between data
        points. The diagonal is zeros.
        
    Returns
    -------
    np.ndarray
        The linkage matrix, as specified in scipy. Briefly, this should be a 2d matrix
        each row containing 4 elements. The first and second element should denote the
        cluster IDs being merged, the third element should be the distance, and the
        fourth element should be the number of elements within this new cluster. Any
        new cluster should be assigned an incrementing ID, e.g. after the first step
        where the first two points are merged into a cluster, the new cluster ID should
        be N, then N+1, N+2, ... in subsequent steps.
    
    """
    import numpy as np
    a=distances
    ban_list=[]
    #print(dd)
    minn=[-1,0,0]
    a_single_ref=np.ones(len(a[0]), dtype=float)
    ret_list=[]
    for _ in range(len(a[0])-1):
        ret_list_helper=[]
        #finding min
        for i in range(len(a)-1):
            if i not in ban_list:
                for j in range(i+1,len(a)):
                    if j not in ban_list:
                        if(minn[0]<0):
                            minn=[a[i][j],i,j]
                        elif(minn[0]>a[i][j]): # stavih = samo da ispadne isto kak u primero posle da se brise
                            minn=[a[i][j],i,j]
    
        ret_list_helper=[minn[1],minn[2],minn[0]]
        #print("minn=",minn)
        ban_list.append(minn[1])
        #print("ban list=",ban_list," i=",i)
        ban_list.append(minn[2])
        #print("ban list=",ban_list," j=",j)
        k=a_single_ref[minn[1]]+a_single_ref[minn[2]]
        a_single_ref=np.append(a_single_ref,k)
        ret_list_helper.append(a_single_ref[len(a_single_ref)-1])
        #print(ret_list_helper)
        ret_list.append(ret_list_helper)
        #print(a_single_ref)
        #print("ban list=",ban_list)
        #print("a=",a)
        aa=np.zeros((len(a[0])+1,len(a[0])+1))
        aa[:len(a),:len(a)]=a
        a=aa

        #print("a=",a)
        #new matrix using a, min and a_ref
        for i in range(0,len(a)-1):
            if i not in ban_list:
                try:
                    a[i][len(a)-1]=(a[i][minn[1]]*a_single_ref[minn[1]]+a[i][minn[2]]*a_single_ref[minn[2]])/(a_single_ref[minn[1]]+a_single_ref[minn[2]])
                except:
                    print("exception")
                    #pass
                try:
                    a[len(a)-1][i]=(a[i][minn[1]]*a_single_ref[minn[1]]+a[i][minn[2]]*a_single_ref[minn[2]])/(a_single_ref[minn[1]]+a_single_ref[minn[2]])
                except:
                    print("exception")
                    #pass
        minn=[-1,0,0]
        #print(a)

    ret_list = np.asarray(ret_list)
    #print(ret_list)
    #<class 'numpy.ndarray'>
    #print(type(ret_list))
    #raise NotImplementedError()
    return ret_list
        
        
def jukes_cantor(p: float) -> float:
    """The Jukes-Cantor correction for estimating genetic distances.
    
    Parameters
    ----------
    p: float
        The proportional distance, i.e. the number of of mismatching symbols (Hamming
        distance) divided by the total sequence length.
        
    Returns
    -------
    float
        The corrected genetic distance.
    
    """
    #import math
    import numpy as np
    return (-3/4*np.log(1 - 4/3*p))
