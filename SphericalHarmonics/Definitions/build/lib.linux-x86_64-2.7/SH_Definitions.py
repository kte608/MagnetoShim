from scipy.special import lpmn
from scipy import sin,cos,pi,sqrt
try:
    from scipy import factorial
except:
    from scipy.misc import factorial
from math import pow,atan2

from KGeometry import Vector
#This should really be replaced with maxima
#import Mathematica
#mat=Mathematica.Mathematica()

def deltanm(n,m):
    """
    The kronecker delta function
    """
    if n==m:
        return 1
    else:
        return 0

def epsilonm(m):
    """
    Neumann Factor from romeo.discretecoildesign.pdf
    """
    return 2-deltanm(m,0)

def factorialQuotient(numerator,denominator):
    """
    result=numerator!/(denominator!)
    """
    diff=numerator-denominator
    if diff>0:
        result=factorial(diff)
    else:
        result=1.0/factorial(-diff)
    return result

def legendrePnm_theta(n,m,theta,harmtype="_{n,m}"):
    """
    _{n}^{m} is the P_{n}^{m} type harmonic and
    _{n,m} is the P_{n,m} type harmonic
    
    """
    if m>n:
        result=0.0
    elif n<0 or m<0:
        # This should really be replaced with maxima
        result=mat.execute_numerical('LegendreP',n,m,cos(theta))
    else:
        result=lpmn(m,n,cos(theta))[0][-1][-1]
    if harmtype=='_{n,m}':
        #print "OldResult: ",result
        result=result*pow(-1.0,float(m))
        #print "newresult: ",result
    return result


def legendrePnm_z(n,m,z):
    if m>n:
        return 0.0
    elif n<0 or m<0:
        # This should really be replaced with maxima
        return mat.execute_numerical('LegendreP',n,m,z)
    else:
        try:
            result=lpmn(m,n,z)[0][-1][-1]
        except ValueError:
            print n," ",m
            raise Exception, "Value Error"
            
        return result 

class SphereHarmonic:
    def __init__(self,n,m,A,B,mode="NONE",harmtype='_{n,m}',Rref=1.0):
        """
        _{n}^{m} is the P_{n}^{m} type harmonic and
        _{n,m} is the P_{n,m} type harmonic

        """
        self.Rref=Rref
        self.n=n
        self.m=m
        self.harmtype=harmtype
        if mode=="NONE":
            print "Invalid type!"
            raise exception
        if mode=="COS_SIN":  #see SHID/Documents/BasisSets/SinCosExp.lyx
            self.cosmag=A
            self.sinmag=B
            self.mag=sqrt(A*A+B*B)
            self.phase=atan2(-B,A)
        if mode=="MAG_PHASE":  #see SHID/Documents/BasisSets/SinCosExp.lyx
            self.mag=A            
            self.phase=B
            self.cosmag=A*cos(B)
            self.sinmag=-A*sin(B)
    def val_point(self,point):
        if self.sinmag!=0.0:
            v1=self.SphericalHarmonicOdd(self.n,self.m,point.rsphere(),point.theta(),point.phi())
        else:
            v1=0.0
        if self.cosmag!=0.0:
            v0=self.SphericalHarmonicEven(self.n,self.m,point.rsphere(),point.theta(),point.phi())
        else:
            v0=0.0
        
        result= self.sinmag*v1+self.cosmag*v0
        #print point.__str__(Mode="SPHERICAL")," ",result
        return result

    def SphericalHarmonicOdd(self,n,m,r,theta,phi):
        #print "n:",n," m:",m," theta:",theta," phi:",phi
        result=pow(r/self.Rref,n)*legendrePnm_theta(n,m,theta,self.harmtype)*sin(m*phi)
        return result

    def SphericalHarmonicEven(self,n,m,r,theta,phi):
        #print "n:",n," m:",m," theta:",theta," phi:",phi
        result=pow(r/self.Rref,n)*legendrePnm_theta(n,m,theta,self.harmtype)*cos(m*phi)
        return result


if __name__ == "__main__":
    # We are going to plot some legendre polynomials
    from numpy import *
    import pylab 
    SH=SphereHarmonic(0,0,0)
    npoints=100
    thetapoints=arange(npoints)*2*pi/npoints
    #for n in range(1,4):
    n=2
    m=2
    y=[SH.legendrePnm(n,m,thetapoint) for thetapoint in thetapoints]       
    pylab.plot(thetapoints,y)
    pylab.show()
