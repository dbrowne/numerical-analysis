# Root finding algorithms for nonlinear equations in one variable
# Author: David Palma
# Ref: Numerical Analysis (Burden and Faires, 9 ed.)
# Last revised: 2011.10.04

import math

# Bisection method for finding roots
# Inputs: function, endpoints, maximum iterations, tolerance
def bisection(f, a, b, n, tol):
    print "Using Bisection Method..."
    print "a = %.9f, b = %.9f" %(a,b)
    i = 1
    FA = f(a)
    while (i <= n):
        p = a + (b - a)/2
        print "Current approximation is %.9f" % p
        FP = f(p)
        if FP == 0 or (b-a)/2 < tol:
            print "Approximate root is %.9f (found in %i steps)" % (p,i)
            return
        i += 1
        if (FA*FP > 0):
            a = p
        else:
            b = p
    print "Root not found in %n steps..."
    print "Stopped at %.9f" % p

# Fixed Point Iteration Method for finding roots
# Inputs: function, initial approximation, maximum iterations, tolerance
def fixed_point(f, p0, n, tol):
    print "Using fixed point iteration..."
    print "p0 = %.9f" % p0
    i = 1
    while (i <= n):
        p = f(p0)
        print "Current approximation is %.9f" % p
        if (abs(p - p0) < tol):
            print "Approximate root is %.9f (found in %i steps)" %(p,i)
            return
        i += 1
        p0 = p
    print "Root not found in %n steps..."
    print "Stopped at %.9f" % p

# Newton's Method for finding roots
# Inputs: function, derivative, initial approximation, maximum iterations, tolerance
def newton_raphson(f, df, p0, n, tol):
    print "Using Newton/Raphson Method..."
    print "p0 = %.9f" % p0
    i = 1
    while (i <= n):
        p = p0 - f(p0)/df(p0)
        print "Current approximation is %.9f" % p
        if (abs(p - p0) < tol):
            print "Approximate root is %.9f (found in %i steps)" %(p,i)
            return
        i += 1
        p0 = p
    print "Root not found in %n steps..."
    print "Stopped at %.9f" % p

# Secant Method for finding roots
# Inputs: function, initial approximations, maximum iterations, tolerance
def secant(f, p0, p1, n, tol):
    print "Using Secant Method..."
    print "p0 = %.9f p1 = %.9f" % (p0,p1)
    q0 = f(p0)
    q1 = f(p1)
    i = 2
    while (i <= n):
        p = p1 - q1*(p1-p0)/(q1-q0)
        print "Current approximation is %.9f" % p
        if (abs(p - p1) < tol):
            print "Root is %.9f (%d iterations)" % (p,i+1)
            return
        i += 1
        p0 = p1
        q0 = q1
        p1 = p
        q1 = f(p)
    print "Root not found in %n steps..."
    print "Stopped at %.9f" % p
	
# Regula Falsi Method for finding roots
# Inputs: function, initial approximations, maximum iterations, tolerance
def regula_falsi(f, p0, p1, n, tol):
    print "Using Regula Falsi (false position) method..."
    print "p0 = %.9f p1 = %.9f" % (p0,p1)
    q0 = f(p0)
    q1 = f(p1)
    i = 2
    while (i <= n):
        p = p1 - q1*(p1-p0)/(q1-q0)
        print "Current approximation is %.9f" % p
        if (abs(p - p1) < tol):
            print "Root is %.9f (%d iterations)" % (p,i+1)
            return
        i += 1
        q = f(p)
        if (q*q1 < 0):
            p0 = p1
            q0 = q1
        p1 = p
        q1 = q
    print "Root not found in %n steps..."
    print "Stopped at %.9f" % p

# Modified Newton's Method for finding roots
# Inputs: function, derivative, second derivative, initial approximation, maximum iterations, tolerance
def modified_newton(f, df, d2f, p0, n, tol):
    print "Using modified Newton's Method..."
    print "p0 = %.9f" %p0

    def g(y): return f(y)/df(y)
    def dg(y): return 1 - (f(y)*d2f(y))/(df(y)**2)

    i = 1
    while (i <= n):
        #p = p0 - (f(p0)*df(p0))/(df(p0)**2 - f(p0)*d2f(p0))
        p = p0 - g(p0)/dg(p0)
        print "Current approximation is %.9f" % p
        if (abs(p - p0) < tol):
            print "Approximate root is %.9f (found in %i steps)" %(p,i)
            return
        i += 1
        p0 = p
    print "Root not found in %n steps..."
    print "Stopped at %.9f" % p

# Aitken's Method for accelerating convergence
# Inputs: function, initial approximation, number of iterations
def aitken(f, p0, n):
    print "Using Aitken's Method..."
    print "p0 = %.9f" % p0
    print ""
    pm = p0
    m = n + 2
    # create array of size m
    arr = [0] * m
    # generate original sequence
    print "Generating original sequence:" 
    i = 1
    while (i <= m):
        arr[i-1] = f(pm)
        q = arr[i-1]
        print "Iteration %d: %.9f" % (i,q)
        pm = q
        i += 1
    # generate accelerated sequence
    print ""
    print "Generating accelerated sequence:"
    pn = p0
    i = 1
    while (i <= n):
        p1 = arr[i-1]
        p2 = arr[i]
        if (p2-2*p1+p0 == 0):
            print "Zero in the denominator."
            print "Stopped at %.9f" % p
            return
        pn = p0 - ((p1-p0)**2)/(p2-2*p1+p0)
        print "Iteration %d: %.9f" % (i,pn)
        p0 = p1
        i += 1

# Steffensen's Method for finding roots
# Inputs: function, initial approximation, maximum iterations, tolerance
def steffensen(f, p0, n, tol):
    print "Using Steffensen's Method..."
    print "p0 = %.9f" % p0
    p = p0
    i = 1
    while (i <= n):
        p1 = f(p0)
        p2 = f(p1)
        if (p2-2*p1+p0 == 0):
            print "Zero in the denominator."
            print "Stopped at %.9f" % p
            return
        p = p0 - ((p1-p0)**2)/(p2-2*p1+p0)
        print "Current approximation is %.9f" % p
        if (abs(p - p0) < tol):
            print "Root is %.9f (%d iterations)" % (p,i)
            return
        i += 1
        p0 = p
    print "Root not found in %n steps..."
    print "Stopped at %.9f" % p

# Horner's Method

# Muller's Method
# Inputs: 
def muller(f, p0, p1, p2, n, tol):
    print "Using Muller's Method..."
    print "p0 = %.9f, p1 = %.9f, p2 = %.9f" % (p0, p1, p2)
    h1 = p1 - p0
    h2 = p2 - p1
    d1 = (f(p1) - f(p0)) / h1
    d2 = (f(p2) - f(p1)) / h2
    d = (d2 - d1)/(h2 + h1)
    i = 3
    while (i <=  n):
        b = d2 + h2*d
        dd = (b*b - 4*f(p2)*d)**0.5
        if (abs(b-dd) < abs(b + dd)):
            e = b + dd
        else:
            e = b - dd
        h = -2*f(p2) / e
        p = p2 + h
        if (abs(h) < tol):
            print "Current approximation is %.9f" % p
            return
        p0 = p1
        p1 = p2
        h1 = p1 - p0
        h2 = p2 - p1
        d1 = (f(p1) - f(p0)) / h1
        d2 = (f(p2) - f(p1)) / h2
        d = (d2 - d1)/(h2 + h1)
        i += 1         
    print "Root not found in %n steps..."
    print "Stopped at %.9f" % p

# Function definitions
def f(x): return 0.5*(x + (3/x))
def df(x): return 2*(math.e)**(-2*x) * ((math.e)**x + 1) * (x*(math.e)**x - 1)
def d2f(x): return 2*(math.e)**(-2*x) * ((x+1)*(math.e)**(3*x) + (math.e)**x + (math.e)**(2*x) + 2)

# Run analysis...
#bisection(f, 1.0, 2.0, 30, 0.0001)
#fixed_point(f, 1.5, 30, 0.0001)
#newton_raphson(f, df, 0.5, 30, 0.00001)
#secant(f, math.e, 4, 30, 0.00001)
#regula_falsi(f, 0.0, 1.0, 30, 0.000001)
#modified_newton(f, df, d2f, 0.5, 30, 0.00001)
#aitken(f, 0.5, 6)
#steffensen(f, 0.0, 30, 0.00001)
#muller(f, 0.0, 0.5, 1.0, 30, 0.00001)
