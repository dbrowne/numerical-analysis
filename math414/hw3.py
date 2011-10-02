import math

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
    print "Root not in %n steps..."
    print "Stopped at %.9f" % p

def secant(f, p0, p1, n, tol):
    print "Finding zeroes using secant method..."
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
    print "Root not in %n steps..."
    print "Stopped at %.9f" % p
	
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
    print "Root not in %n steps..."
    print "Stopped at %.9f" % p

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
    print "Root not in %n steps..."
    print "Stopped at %.9f" % p

def f(x): return x**2 - 2*x*(math.e)**(-x) + (math.e)**(-2*x) 
#def f(x): return math.cos(x + math.sqrt(2)) + x*((x/2) + math.sqrt(2))
#def df(x): return 2*x - 2*math.e**(-2*x) - 2*math.e**(-x) + 2*x*math.e**x
def df(x): return 2*(math.e)**(-2*x) * ((math.e)**x + 1) * (x*(math.e)**x - 1)
#def df(x): return x - math.sin(x + math.sqrt(2)) + math.sqrt(2)
#def d2f(x): return 2 + 4*math.e**(-2*x) + 2*math.e**(-x) + 2*x*math.e**x + 2*math.e**x
def d2f(x): return 2*(math.e)**(-2*x) * ((x+1)*(math.e)**(3*x) + (math.e)**x + (math.e)**(2*x) + 2)

newton_raphson(f, df, 0.5, 30, 0.00001)
modified_newton(f, df, d2f, 0.5, 30, 0.00001)
#secant(f, math.e, 4, 30, 0.00001)
#regula_falsi(f, 0.0, 1.0, 30, 0.000001)
