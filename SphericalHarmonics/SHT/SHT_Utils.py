#!/usr/bin/env python
from SH_Definitions import SphereHarmonic
from SHT_Points import *



def harmonizer(meas_points,
               harm_coeffs,
               mode="NONE",
               pointType="THETAPHI",
               harmtype='_{n}^{m}'):
    """This function takes a number of measurement points (meas_points)
    as a list of the form ((theta,phi),(theta,phi),...) that describes
    the points where the function is to be calculated. It also takes
    a list of harmonic coefficients in one of two forms:
    'COS_SIN':
    ((n,m,cosmag,sinmag),(n,m,cosmag,sinmag),etc)
    'MAG_PHASE':
    ((n,m,mag,phase),(n,m,mag,phase),etc)
    of the function to be calculated at those points. The result is a list of
    two lists. The first is of the form:
    ((theta,phi,realvalue),...)
    The second is of the form (realvalue,realvalue,...)
    There are no imaginary values computed 
    """
    result=[]
    rdata=[]
    for mp in meas_points:
        if pointType=="THETAPHI":
            theta=mp[0]
            phi=mp[1]
            p=Vector((1.0,theta,phi),type="SPHERICAL")
        elif pointType=="VECTOR":
            p=mp
            theta=p.theta()
            phi=p.phi()
        sumr=0.0
        for (n,m),(A,B) in harm_coeffs:
            
            #print type(n)
            #coeff=harmcoeff(n,m)

            if mode=="NONE":
                print "Must Specify type"
                raise Exception
            if mode=="COS_SIN":
                cosmag=A
                sinmag=B
                val=SphereHarmonic(n,m,cosmag,sinmag,
                                   mode="COS_SIN",
                                   harmtype=harmtype).val_point(p)
            if mode=="MAG_PHASE":
                mag=A
                phase=B
                val=SphereHarmonic(n,m,mag,phase,mode="MAG_PHASE",harmtype=harmtype).val_point(p)

            sumr=sumr+val
        rdata.append(sumr)
        result.append((theta,phi,sumr))
    return result,rdata



def spherefieldplotterMAT(bw,fieldvalues,type="2D",title=''):
    """
    I assume that fieldvalues is a list [val,val,val,...] where
    each value corresponds to a field point in the list:
    fieldpointsVect=[Vector((Rsphere,theta,phi),type="SPHERICAL") for theta,phi in measurement_points(bw)]
    """
    import pylab as p
    

    fieldpointsVect=[Vector((1.0,theta,phi),type="SPHERICAL") for theta,phi in measurement_points(bw)]
    
    x=azimuthal_measurement_points(bw)


    y=declination_measurement_points(bw)

 
    
    Z=[]
    i=0
    for theta in y:
        tmp=[]
        for phi in x:
            if fieldpointsVect[i].theta() != theta:
                print "i: ",i," Bad Theta"
                print theta
                print fieldpointsVect[i].theta()
                print phi
                print fieldpointsVect[i].phi()
                for p in fieldpointsVect:
                    print "theta: ",p.theta()
                    print "phi:   ",p.phi()
                print "thetas: ",x
                print "phis: ",y
            
                raise Exception
            if fieldpointsVect[i].phi() != phi:
                print "i: ",i," Bad Phi"
                print theta
                print fieldpointsVect[i].theta()
                print phi
                print fieldpointsVect[i].phi()
                raise Exception
            tmp.append(fieldvalues[i])
            i=i+1
        tmp.append(tmp[0])
        Z.append(array(tmp))
    Z.append(array(Z[0]))
    Z=array(Z)
    #print len(Z)
    #print len(Z[0])


    y.append(pi)
    x.append(2*pi)
    y=array(y)
    x=array(x)
    #print "xlen: ",len(x),"ylen: ",len(y)
    X,Y=p.meshgrid(x,y)

    if type=="2D":
        p.hold(True)
        p.imshow(Z,origin='lower')
        p.contour(Z,10,origin='lower')
        p.xlabel('phi (azimuth)')
        p.ylabel('theta (declination)')
        p.title(title)
        p.show()
    if type=="3D":
        import matplotlib.axes3d as p3
        fig=p.figure()
        ax=p3.Axes3D(fig)
        ax.contour3D(X,Y,Z,50)
        ax.set_xlabel('phi (azimuth)')
        ax.set_ylabel('theta (declination)')
    
        #fig.add_axes(ax)
        p.show()

def harmplotterMAT(bw,harms,Rsphere,mode="COS_SIN",type="2D",title='',harmtype='_{m}^{m}'):
    """Matplot some harmonics"""
    import pylab as p
    
    x=azimuthal_measurement_points(bw)
    x.append(2*pi)
    x=array(x)
    y=declination_measurement_points(bw)
    y.append(pi)
    y=array(y)
    X,Y=p.meshgrid(x,y)
    Z=[]
    for Xa,Ya in zip(X,Y):
        tmp=[]
        for xa,ya in zip(Xa,Ya):
            fp=Vector([Rsphere,ya,xa],type="SPHERICAL")
            d,rdata=harmonizer([fp],harms,mode=mode,pointType="VECTOR",harmtype=harmtype)
            result=float(rdata[0])
            tmp.append(result)
        Z.append(array(tmp))
    Z=array(Z)
    print len(Z)
    print len(Z[0])
    if type=="3DSPHERE":
        import matplotlib.axes3d as p3
        print "I am not sure if this is working properly"
        print "See scipy wiki cookbook for mplot3D"
        from pprint import pprint
        phis=y#r_[0:2*pi:100j]
        thetas=x#r_[0:pi:100j]
    
    
        xvect=1*outer(cos(phis),sin(thetas))
        yvect=1*outer(sin(phis),sin(thetas))
        zvect=1*outer(ones(size(phis)),cos(thetas))
        for i,j in zip(range(len(thetas)),range(len(phis))):
            print Z[i][j]
            xvect[i][j]=xvect[i][j]*Z[i][j]
            yvect[i][j]=yvect[i][j]*Z[i][j]
            zvect[i][j]=zvect[i][j]*Z[i][j]
        #print(len(xvect))
        #print(len(xvect[0][0]))
        fig=p.figure()
        ax=p3.Axes3D(fig)
        ax.plot_surface(xvect,yvect,zvect)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        p.show()
    if type=="3D":
        fig=p.figure()
        ax=p3.Axes3D(fig)
        ax.contour3D(X,Y,Z,50)
        ax.set_xlabel('phi (azimuth)')
        ax.set_ylabel('theta (declination)')
    
        #fig.add_axes(ax)
        p.show()
    if type=="2D":
        p.hold(True)
        p.imshow(Z,origin='lower')#,origin='lower',extent=[-1,1,-1,1])
   

        p.contour(Z,10,origin='lower')#y,x,Z,10)#,origin='lower',extent=[-1,1,-1,1])
        p.xlabel('phi (azimuth)')
        p.ylabel('theta (declination)')
        p.title(title)
        p.show()

if __name__ == "__main__":
    bw=6
    Rsphere=1.0
    fieldpointsVect=[Vector((Rsphere,theta,phi),type="SPHERICAL") for theta,phi in measurement_points(bw)]
    fieldvalues=[4.0 for f in fieldpointsVect]
    harmplotterMAT(bw,[((1,1),(1,0))],Rsphere)
    spherefieldplotterMAT(bw,fieldvalues,type="2D",title='')
