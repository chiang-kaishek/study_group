__doc__ = r"""
Optimiaztion linear and quadratic cost, EXPLOIT KKT STRUCTURE, no trade bound position bound
Author: ChiKong
Platform: python 2.7 with numpy, scipy and cvxopt
Version: V1.0
Date 2017-03-06
"""

import numpy
from cvxopt import matrix, solvers, spmatrix,lapack,blas, uniform, mul, spdiag
import time
#import opt_lqtc_1p

solvers.options['show_progress'] = False

def opt_lqtc_1p_exkkt( S , r  , w0 = None , c = None  , D = None  ,    **kwargs): 
    """
    Solve the optimization problem 

    minimize     1/2 w'*S*w - r' * w + c' * |w-w0| + (w-w0)'D(w-w0)
    subject to   no constraint 
    
    Input Argument( basic usage ).

                 S : covariance 
                 r : expected return

    Optional Argument.
                 w0: Initial Position
                 c : Linear tcost constant
                 D : Quadratic tcost constant

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
    #if w0 is None: 
    #    w0     =  matrix( 0.0 , (dim,1))
    #if c is None: 
    #    c      =  matrix( 0.0 , (dim,1))
    #if D is None: 
    #    D      =  matrix( 0.0 , (dim,dim))
    w0 = matrix(w0) if w0 is not None else matrix(0., (dim, 1))
    c = matrix(c) if c is not None else matrix(0., (dim, 1))
    D = matrix(D) if D is not None else matrix(0., (dim, dim))
            
    #  formulation of bigger optimization problem  
    P     =   matrix([[ S + 2*D , -S + 2*D],[-S + 2*D , S + 2*D]])
    q     =   matrix( [S*w0 - r + c  , -S*w0 + r + c  ] )

    #S1 = S + 2 * D
    #S2 = -S + 2 * D
    
    def G(x,y,alpha = 1.0 , beta = 0.0, trans = 'N'):
    # funciton valued constraint
        if trans == 'N':
            # y:= alpha * G * x + beta * y
            y[:]    =  - alpha * x +  beta * y   
        else:
            # y:= alpha * G' * x + beta*y
            y[:]    =  - alpha * x +  beta * y  
            
    def Px(x,y,alpha = 1.0,beta=0.0):
        # function valued matrix
        y[:]   =   alpha * P * x + beta * y
        

    def F(W):
        """
        Return a function f(x,y,z) that solves
        [P , G' W^-1]  [ux]       [bx]
        [G ,  -W    ]  [uy]    =  [bz]

        """
        #d   =    spdiag(matrix(numpy.array(W['d'])))
        #dinv=    spdiag(matrix(numpy.array(W['di'])))
        d   =    spdiag(W['d'])
        dinv=    spdiag(W['di'])

        #KKT1 =    d*( P * d + dinv ) 
        KKT1 =    d*P*d + spdiag(matrix(1.0,(2*dim,1)))
        lapack.potrf(KKT1)
        
        #raw_input('inputpppp')
        def f(x,y,z):
            uz   =  -d* (x + P * z)
            #uz  =   matrix(numpy.linalg.solve(KKT1, uz))  # slow version
            #lapack.gesv(KKT1,uz)  #  JZ: gesv have cond issue 
            lapack.potrs(KKT1,uz)
            x[:]  =   matrix( -z - d *  uz)
            blas.copy(uz , z)

        return f


    h     =   matrix(numpy.bmat([[r*0],[r*0]]))
    dims  =  {'l': 2*dim, 'q': [], 's': []}
    #  solver
    sol   =   solvers.coneqp(Px,q , G,h , dims , kktsolver = F) # exploit structured kkt
    #sol   =   solvers.coneqp(P,q , G,h , dims)
    x     =   sol['x']
    status=   sol['status']
    print status
    return {'x':x, 'b':x[0:dim] ,'s':x[dim:2*dim],'dw':x[0:dim] -x[dim:2*dim],'status':status }








#if __name__ == "__main__":
    #  a running example
    

