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
    print "Root not found to tolerance %.9f in %n..." % tol
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
    print "Root not found to tolerance %.9f in %n..." % tol
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
    print "Root not found to tolerance %.9f in %n..." % tol
    print "Stopped at %.9f" % p
    
#def f(x): return 230*x**4 + 18*x**3 + 9*x**2 - 221*x - 9
def f(x): return (x-2)**2 - math.log(x)
#def df(x): return 920*x**3 + 54*x**2 + 18*x - 221
def df(x): return 2*(x-2)-(1/x)

newton_raphson(f, df, 3.0, 30, 0.00001)
secant(f, math.e, 4, 30, 0.00001)
#regula_falsi(f, 0.0, 1.0, 30, 0.000001)
