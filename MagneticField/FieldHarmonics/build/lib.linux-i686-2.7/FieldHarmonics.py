#!/usr/bin/env python
from SH_Definitions import SphereHarmonic
from KGeometry import Vector

def delta(a,b):
    if a==b:
        return 1
    else:
        return 0

def ScalarPot2BzField(n_scpot,m_scpot,A_scpot,B_scpot):
    
    n=n_scpot-1
    m=m_scpot
    A=A_scpot*(n_scpot+m_scpot)
    B=B_scpot*(n_scpot+m_scpot)
    
    if n_scpot==0:
        n=m=A=B=None
    return (n,m),(A,B)

class FieldHarmonic:
    """
    A Field harmonic described according to Bz=A T_{n,m}+B T_{n,m}^{\prime}
    """
    def __init__(self,n,m,A,B):
        self.n=n
        self.m=m
        self.A=A
        self.B=B
        self.Bz=SphereHarmonic(n,m,A,B,mode="COS_SIN",harmtype='_{n,m}')




        fact=-(1+delta(m,0))/(2.0*(n+m+1))        
        self.Bxa=SphereHarmonic(n,m+1,A*fact,B*fact,mode="COS_SIN",harmtype='_{n,m}')
        self.Bya=SphereHarmonic(n,m+1,-B*fact,A*fact,mode="COS_SIN",harmtype='_{n,m}')
        fact=((1-delta(m,0))*(n+m)*(n+m+1))
        if fact != 0.0:
            fact=fact/(2.0*(n+m+1))
        self.Bxb=SphereHarmonic(n,m-1,fact*A,fact*B,mode="COS_SIN",harmtype='_{n,m}')
        self.Byb=SphereHarmonic(n,m-1,fact*B,-fact*A,mode="COS_SIN",harmtype='_{n,m}')
        
        
    def A_val_point(self,point):
        # The vector potential at the point
        pass

    def S_val_point(self,point):
        # The scalar potential at the point
        pass

    def B_val_point(self,point):
        Bx=self.Bxa.val_point(point)+self.Bxb.val_point(point)
        By=self.Bya.val_point(point)+self.Byb.val_point(point)
        Bz=self.Bz.val_point(point)
        result=Vector((Bx,By,Bz))
        return result

class FieldHarmonicEnsemble:
    def __init__(self,harmonics):
        """
        Harmonics is either a list or a dictionary
        of harmonics (dictionary for now) of the 
        form {(n,m):(Ba,Bb),...} and we can evaluate
        the whole field.
        """
        self.FieldHarmonics=[]
        for (n,m),(Ba,Bb) in harmonics.iteritems():
            self.FieldHarmonics.append(FieldHarmonic(n,m,Ba,Bb))
    def B_val_point(self,point):
        B=Vector((0,0,0))
        for f in self.FieldHarmonics:
            B=B+f.B_val_point(point)
        return B
        

if __name__ == "__main__":
    fh=FieldHarmonic(2,3,1,0)
    print fh.B_val_point(Vector((1,1,1)))
