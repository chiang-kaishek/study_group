__doc__ = r"""
Optimiaztion linear and quadratic cost, slow, non-optimal version, direct implementation
Author: ChiKong
Platform: python 2.7 with numpy, scipy and cvxopt
Version: V1.0
Date 2017-03-06
"""


import numpy
from cvxopt import matrix, solvers, spmatrix
import time

def opt_lqtc_1p( S , r  , w0 = None , c = None  , D = None , pl = None , pu = None , tl = None , tu = None ,    **kwargs): 
    """
    Solve the optimization problem 

    minimize     1/2 w'*S*w - r' * w + c' * |w-w0| + (w-w0)'D(w-w0)
    subject to   pl  <=  w  <=   pu
                 tu <= (w-w0) <= tl
    
    Input Argument( basic usage ).

                 S : covariance 
                 r : expected return

    Optional Argument.
                 w0: Initial Position
                 c : Linear tcost constant
                 D : Quadratic tcost constant
                 tl,tu:  lower and upper trade bound
                 pl,pu:  lower and upper position bound

    Output Argument
                 x : solution to qp
                 b : optimal buys
                 s : optimal sells
                 dw: optimal trades
    """ 
      
    S         =     matrix(S)
    r         =     matrix(r)
    dim       =     S.size[0]
    max_bound =  1e6
    if w0 is None: 
        w0     =  matrix( 0.0 , (dim,1))
    if c is None: 
        c      =  matrix( 0.0 , (dim,1))
    if D is None: 
        D      =  matrix( 0.0 , (dim,dim))
    if pl is None: 
        pl     =  matrix( -max_bound, (dim,1))
    if pu is None: 
        pu     =  matrix( max_bound , (dim,1))
    if tl is None: 
        tl     =  matrix( -max_bound, (dim,1))
    if tu is None: 
        tu     =  matrix( max_bound, (dim,1))

    #  formulation of bigger optimization problem  
    P     =   matrix([[ S + 2*D , -S + 2*D],[-S + 2*D , S + 2*D]])
    q     =   matrix( [S*w0 - r + c  , -S*w0 + r + c  ] )
    I     =   spmatrix(1.0,range(dim),range(dim))
    G     =   matrix(  [[I , -I , -I , 0*I ] ,[-I , I , 0*I , -I]]) 
    h1    =   numpy.minimum(pu - w0 , tu)
    h2    =   numpy.minimum(-(pl - w0) , -tl)
    h     =   matrix(numpy.bmat([[h1],[h2],[r*0],[r*0]]))
    #  solver
    sol   =   solvers.qp(P,q , G,h)
    x     =   sol['x']
    status=   sol['status']
    print status
    return {'x':x, 'b':x[0:dim] ,'s':x[dim:2*dim],'dw':x[0:dim] -x[dim:2*dim],'status':status }




if __name__ == "__main__":
    #  a running example
    