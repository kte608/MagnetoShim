#!/usr/bin/env python
import bSHT
from math import pi,sqrt,atan2,pow
from KStandard import inclusive_range
#from SH_Definitions import *
from KGeometry import *
from SHT_Points import *
def SHT(bw,
        realdata=None,
        imagdata=None,
        mode="NONE",
        normalization="NONE",
        Rref=None,
        harmtype='_{n,m}'):
    """
    if the normalization == "NONE", as is the default, then the harmonics are simply multiplied by unity

    if Rref is given the coefficients will be divided by Rref^n so that they will not depend on the
    radius of the sphere on which they were taken.

   _{n}^{m} is the P_{n}^{m} type harmonic and
   _{n,m} is the P_{n,m} type harmonic

    if mode == "NORMAL":
        # Real and imaginary data are given
        #
        #
        # It is important to realize the relationship (given in s2kit_fx.pdf) between coefficients where:
        # f(l,m)=(-1)^m Conjugate(f(l,m))

        
    if mode == "COS_SIN":
        # We only take real data as an input and we return a list
        # [((n,m),(cosinemag,sinemag)),...]
        #
        # where the realdata fits FUNC=SUM[n,m,Normalization*sphereharmonic*(cosinemag*cos()+sinemag*sin())]        
        #

    if mode == "COS_SIN_NOWRAP":
        # we only take real data as an input and we return a list
        # as in COS_SIN but it does not do any fancy wrapping of the negative harmonics

    if mode == "MAG_PHASE":
        # We only take real data as an input and we return a list
        # [((n,m),(mag,phase)),....]
        #
        # where the realdata fits FUNC=SUM[n,m,Normalization*sphereharmonic*(mag*cos(m*(phi-psi)))]        
        #

        """

    if Rref == None:
        print "Are you sure you don't want to use Rref??!!??"
        raise Exception

    if mode == "NONE":
        print "You have called SHT incorrectly."
        return
    if mode == "NORMAL":
        # Real and imaginary data are given
        #
        #
        # It is important to realize the relationship (given in s2kit_fx.pdf) between coefficients where:
        # f(l,m)=(-1)^m Conjugate(f(l,m))
        if Rref != None:
            print "Rref is not considered in the Normal Mode yet"
            raise Exception
        return SHT_ri(bw,realdata,imagdata,normalization=normalization,harmtype=harmtype)
    if mode == "COS_SIN":
        # We only take real data as an input and we return a list
        # [((n,m),(cosinemag,sinemag)),...]
        #
        # where the realdata fits FUNC=SUM[n,m,Normalization*sphereharmonic*(cosinemag*cos()+sinemag*sin())]        
        #
        if imagdata != None:
            print "You cannot imput imaginary data in this mode: COSINE_SINE"
            return
        imagdata=[0 for i in realdata]
        resultlist=SHT_ri(bw,realdata,imagdata,normalization=normalization,harmtype=harmtype)
        # For the zonal harmonics only the cosinemagnitude makes any sense and it equals the magnitude
        # of the normal e^iphi version of the harmonic transform.
        #

        # 
        result=[]
        
        for i in resultlist[:len(genindicies_top(bw))]:
            m=i[0][1]
            if m==0:
                #Zonal harmonics
                if Rref != None:
                    n=i[0][0]
                    magcos=i[1][0]
                    magsin=i[1][1]
                    magcos=magcos/pow(Rref,n)
                    magsin=magsin/pow(Rref,n)
                result.append((i[0],(magcos,magsin)))
            else:
                #Tesseral harmonics need some work. See SinCosExp.lyx under SHID/Documents/BasisSets
                n=i[0][0]
                real=i[1][0]
                imag=i[1][1]
                a=sqrt(real*real+imag*imag)
                A=atan2(imag,real) 
                cosmag=2*a*cos(A)
                sinmag=-2*a*sin(A)
                if Rref != None:
                    #print "I am dividing by Rref"
                    cosmag=cosmag/pow(Rref,n)
                    sinmag=sinmag/pow(Rref,n)
                result.append((i[0],(cosmag,sinmag)))

        
      
        return result
    if mode == "MAG_PHASE":
        # We only take real data as an input and we return a list
        # [((n,m),(mag,phase)),....]
        #
        # where the realdata fits FUNC=SUM[n,m,Normalization*sphereharmonic*(mag*cos(m*(phi-psi)))]        
        #
        if Rref != None:
            print "Rref is not considered in the MAG_PHASE Mode yet"
            raise Exception
        if imagdata != None:
            print "You cannot imput imaginary data in this mode"
            return
        imagdata=[0 for i in realdata]
        resultlist=SHT_ri(bw,realdata,imagdata,normalization=normalization)
        result=[]
        for i in resultlist[:len(genindicies_top(bw))]:
            m=i[0][1]
            if m==0:
                #Zonal harmonics need no alteration
                result.append(i)
            else:
                #Tesseral harmonics need some work. See SinCosExp.lyx under SHID/Documents/BasisSets
                real=i[1][0]
                imag=i[1][1]
                a=sqrt(real*real+imag*imag)
                A=atan2(imag,real) 
                mag=2*a
                phase=A
                result.append((i[0],(mag,phase)))
        return result



def SHT_ri(bw,realdata,imagdata,normalization="NONE",harmtype='_{n}^{m}'):
    """Given the bandwidth of the problem (the order that is one larger
    than the maximum order considered) and real and imaginary data
    lists measured or calculated on the points defined by measurement_points
    then this function returns a list:
    [((n,m),(realcoeff,imagcoeff)),etc]
    which describes the coefficients of the spherical harmonics found
    in the data.

   _{n}^{m} is the P_{n}^{m} type harmonic and
   _{n,m} is the P_{n,m} type harmonic

    
    """
    #print "I am in the python part of the SHT"
    if normalization != "NONE":
        print "I have not dealt with that situation yet"
        return
    numdecl=len(declination_measurement_points(bw))
    numaz=len(azimuthal_measurement_points(bw))
    num=numaz*numdecl
    if(len(realdata)!=num or len(imagdata)!=num): return None
    tmpresult= zip(genindicies(bw),bSHT.SHT((bw,realdata,imagdata)))
    result=[]
    if harmtype=='_{n,m}':
        for (n,m),(A,B) in tmpresult:
            result.append([(n,m),(pow(-1,m)*A,pow(-1,m)*B)])
    else:
        result=tmpresult
    return result
    

def factorial(n):
    """Computes the factorial of n"""
    if n==0: return 1

    result=1
    for i in inclusive_range(2,n):
        result=result*i
    return result

def harmcoeff(n,m):
    """Computes the normalization coefficient for harmonic n,m"""
    result = sqrt(((2*n+1)*factorial(n-m))/((4*pi)*factorial(n+m)))
    return result




        
if __name__ == "__main__":
    print declination_measurement_points(3)
    print azimuthal_measurement_points(3)
    
