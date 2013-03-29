#!/usr/bin/env python
# This is a code for geometry.
# author: Karl Edler

##############
# TODO
##############
#
# Find out about 4x4 matricies and how they can be used for
# rotations,translations,scales,and skews.
#
# Might want to integrate what is here into what has been written in
# ScientificPython Geometry (not scipy)
#
# Something must be done to address the difference between a coordinate system
# change and an object change (rotation or translation)
#

from math import *
from KStandard import BaseClass,basicallyzero,attributesFromDict
from numpy.linalg import inv as invert
from numpy import array
class Vector(BaseClass):
    """ Documentation for the Vector


    # WARNING WARNING WARNING. There is a conceptual difference between 'Point'
    # and 'Vector' which is not properly addressed in this code.
    # I can define the spherical coordinates of a point but if I define the
    # spherical components of a vector its Cartesian components will depend on
    # position.

    More details
    SPHERICAL:
      r,theta,phi (theta is declination and phi is azimuth) 
    CYLINDRICAL:
      r,phi,z     

    """
    validtypes=set(["CARTESIAN","SPHERICAL","CYLINDRICAL","TOROIDAL"])


    def __init__(self,coords,vtype="CARTESIAN",_ClassName_="vector"):

        if not vtype in Vector.validtypes:
            raise Exception
        self._ClassName_=_ClassName_
        self._rsphere = None
        self._theta = None
        self._phi = None
        self._rcylind = None
        self._magnitude2 = None
        self._uvect = None
        # See Korn and Korn section 6.5-1 for lots
        # of different curvilinear coordinate systems.
        if vtype == "TOROIDAL":
            # See Korn and Korn page 186
            #
            # c and nu set the shape of the toroid.
            # The center of the toroid's ring cross section
            # is offset from the z-axis through the center of
            # the toroid by an amount cT*coth(nuT) and the 
            # radius of this circular cross-section is cT/sinh(nuT)
            #
            # also see KGeometry_ToroidalCoordsToroidShape.jpg in this directory.
            #
            # phi is the angle about the z-axis through the center of the toroid
            # and goes from 0 to 2*pi
            # theta is an angle that goes around the surface of the toroid and
            # also goes from 0 to 2*pi
            #
            self._nuT=coords[0]     # tau   in Korn and Korn
            self._thetaT=coords[1]  # sigma in Korn and Korn
            self._phiT=coords[2]    # phi   in Korn and Korn
            if len(coords)==4:
                self._cT=coords[3]  # a     in Korn and Korn
            else:
                self._cT=1.0
            denom=float(cosh(self._nuT)-cos(self._thetaT))
            self._x=(self._cT*sinh(self._nuT)*cos(self._phiT))/denom
            self._y=(self._cT*sinh(self._nuT)*sin(self._phiT))/denom
            self._z=(self._cT*sin(self._thetaT))/denom

        if vtype == "SPHERICAL":
            self._rsphere=coords[0]
            self._theta=coords[1]
            self._phi=coords[2]
            self._x=self._rsphere*sin(self._theta)*cos(self._phi)
            self._y=self._rsphere*sin(self._theta)*sin(self._phi)
            self._z=self._rsphere*cos(self._theta)
            #
        if vtype == "CARTESIAN":
            self._x=coords[0]
            self._y=coords[1]
            self._z=coords[2]
        if vtype == "CYLINDRICAL":
            self._rcylind=coords[0]
            self._phi=coords[1]
            self._z=coords[2]
            self._x=self._rcylind*cos(self._phi)
            self._y=self._rcylind*sin(self._phi)


    def ZCompare(a,b):
        if a.z() > b.z():
            return 1
        if a.z() < b.z():
            return -1
        else:
            return 0

    def RCylindCompare(a,b):
        if a.rcylind() > b.rcylind():
            return 1
        if a.rcylind() < b.rcylind():
            return -1
        else:
            return 0

    def PhiCompare(a,b):
        if a.phi() > b.phi():
            return 1
        if a.phi() < b.phi():
            return -1
        else:
            return 0

    def magnitude2(self):
        if self._magnitude2 == None:
            self._magnitude2 = self._x*self._x+self._y*self._y+self._z*self._z
        return self._magnitude2

    def magnitude(self):
        return self.rsphere()

    def length(self):
        return self.magnitude()

    def unitvector(self):
        if self._uvect == None:
            l=self.rsphere()
            if l==0.0:
                self._uvect=None
            else:
                self._uvect=Vector((self._x/l,self._y/l,self._z/l))
        return self._uvect
        


    #CARTESIAN
    def x(self):
        return float(self._x)
    def y(self):
        return float(self._y)
    def z(self):
        return float(self._z)
    def xyz(self):
        return (self.x(),self.y(),self.z())
    #SPHERICAL
    def rsphere(self):
        if self._rsphere == None:
            self._rsphere = sqrt(self.magnitude2())
        return self._rsphere

    def theta(self):
        if self._theta == None:
            if self.rsphere()==0.0:
                self._theta=0.0
            else:
                self._theta=acos(self._z/self._rsphere)
        return self._theta

    def phi(self,format="BAL"):
        if self._phi==None:
            self._phi=atan2(self._y,self._x)
        if format=="BAL":
            return self._phi
        if format=="POS":
            phi=self._phi
            if phi<0:
                phi=phi+2*pi
            if phi>2*pi:
                phi=phi-2*pi
            return phi

    #CYLINDRICAL
    def rcylind(self):
        if self._rcylind == None:
            self._rcylind=sqrt(self._x*self._x+self._y*self._y)
        return self._rcylind






    def __str__(self,Mode="CARTESIAN"):
        if Mode == "CARTESIAN":
            return "XYZ:("+self.x().__str__()+","+self.y().__str__()+","+self.z().__str__()+")"
        elif Mode == "CYLINDRICAL":
            return "RPZ:("+self.rcylind().__str__()+","+self.phi().__str__()+","+self.z().__str__()+")"
        elif Mode == "SPHERICAL":
            return "RTP:("+self.rsphere().__str__()+","+self.theta().__str__()+","+self.phi().__str__()+")"
        elif Mode == "TOROIDAL":
            return 
        else:
            print "Unknown Mode"
            raise Exception

    def negate(self):
        return Vector((-self._x,-self._y,-self._z))

    def MirrorYZ(self):
        # Only negate X
        return Vector((-self._x,self._y,self._z))

    def MirrorXZ(self):
        # Only negate Y
        return Vector((self._x,-self._y,self._z))
 
    def MirrorXY(self):
        # Only negate Z
        return Vector((self._x,self._y,-self._z))


    def add(self,othervector):
        x=self._x+othervector._x
        y=self._y+othervector._y
        z=self._z+othervector._z
        result=Vector((x,y,z))
        return result

    def __add__(self,other):
        return self.add(other)

    def __sub__(self,other):
        return self.subtract(other)

    def scale(self,factor):
        x=self._x*factor
        y=self._y*factor
        z=self._z*factor
        return Vector((x,y,z))

    def subtract(self,othervector):
        x=self._x-othervector._x
        y=self._y-othervector._y
        z=self._z-othervector._z
        result=Vector((x,y,z))
        return result

    def dot(self,othervector):
        x=self._x*othervector._x
        y=self._y*othervector._y
        z=self._z*othervector._z
        result=x+y+z
        return result

    def __mul__(self,othervector):
        """
        '*' means dot product
        """
        return self.dot(othervector)

    def cross(self,othervector):
        if othervector==None:
            raise Exception
        rx=self.y()*othervector.z()-self.z()*othervector.y()
        ry=self.z()*othervector.x()-self.x()*othervector.z()
        rz=self.x()*othervector.y()-self.y()*othervector.x()
        result=Vector((rx,ry,rz))
        return result

    def __mod__(self,othervector):
        """
        '%' means cross product
        """
        return self.cross(othervector)
    def translate(self,vect):
        return Vector((self._x+vect._x,
                       self._y+vect._y,
                       self._z+vect._z))      

    def rotate(self,rotmatrix):
        """ rotate the vector by using the transformation matrix"""
        a=rotmatrix
        import scipy
        x_y_z=scipy.dot(rotmatrix,scipy.array([self._x,self._y,self._z]))
        return Vector(x_y_z)


    def transform(self,orient):
        """orient has the form [RotationMatrix,TranslationVector] or
        [TranslationVector,RotationMatrix] depending on the order
        in the list the order of operations is switched"""
        if hasattr(orient[0],'__len__'):
            # The first one is the rotation matrix
            return self.transformRT(orient[0],orient[1])
        else:
            # The first one is the translation vector
            return self.transformTR(orient[0],orient[1])

    def transformRT(self,rotmatrix,transvect):
        """ First Rotate and then translate """
        tmp=self.rotate(rotmatrix)
        result=tmp.translate(transvect)
        return result

    def transformTR(self,transvect,rotmatrix):
        """ First Translate and then Rotate """
        tmp=self.translate(transvect)
        result=tmp.rotate(rotmatrix)
        return result


def find_yaw_pitch_roll_transform(yaw,pitch,roll):
    import scipy
    yt=find_yaw_transform(yaw)
    pt=find_pitch_transform(pitch)
    rt=find_roll_transform(roll)
    tmp=scipy.dot(pt,rt)
    return scipy.dot(yt,tmp)
        
def find_yaw_transform(yaw):
    """ Yaw about the z-axis, yaw is in radians"""
    return [[cos(yaw),sin(yaw),0],
            [-sin(yaw),cos(yaw),0],
            [0,0,1]]

def find_pitch_transform(pitch):
    """ Pitch about the y-axis"""
    return [[cos(pitch),0,-sin(pitch)],
            [0,1,0],
            [sin(pitch),0,cos(pitch)]]

def find_roll_transform(roll):
    """ Roll about the x-axis """
    return [[1,0,0],
            [0,cos(roll),sin(roll)],
            [0,-sin(roll),cos(roll)]]


def find_rotation_about_axis(rotaxis,angle):
    """
    returns the transformation matrix which causes
    a rotation of 'angle' radians about the 'rotaxis'


    see http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
    
    The algorithm (you should simplify it to a single step)
    is:
    (1) Translate space so that the rotation axis passes through the origin.

    (2) Rotate space about the z axis so that the rotation axis lies in the xz plane.

    (3) Rotate space about the y axis so that the rotation axis lies along the z axis.

    (4) Perform the desired rotation by theta about the z axis.

    (5) Apply the inverse of step (3).

    (6) Apply the inverse of step (2).

    (7) Apply the inverse of step (1). 
    """

    pass

def find_rotation_transform(eulerangles):
    """ returns the transformation matrix which rotates the vector by
    (alpha,beta,gamma) which are rotations
    first by alpha about z, then by beta about y', and finally by gamma
    about z'. This is described by Sakurai and in rotations.lyx """
    alpha=eulerangles[0]
    beta=eulerangles[1]
    gamma=eulerangles[2]
    ca=cos(alpha)
    cb=cos(beta)
    cg=cos(gamma)
    sa=sin(alpha)
    sb=sin(beta)
    sg=sin(gamma)

    a00=ca*cb*cg-sa*sg
    a01=ca*cb*sg+cg*sa
    a02=-ca*sb
    a10=-ca*sg-cg*cb*sa
    a11=cg*ca-cb*sg*sa
    a12=sa*sb
    a20=cg*sb
    a21=sg*sb
    a22=cb
    
    return [[a00,a01,a02],
            [a10,a11,a12],
            [a20,a21,a22]]

def find_vector_z_ab(vect):
    """
    find the rotation angles alpha beta determined by assuming
    that a vector (0,0,z) will be rotated to point in the 
    direction of vect. See rotations.lyx for the mathematics.
    """
    xp=vect.x()
    yp=vect.y()
    zp=vect.z()
    x=0
    y=0
    z=1.0
    #beta=acos(zp/float(z))
    alpha=atan2(-yp,float(xp))
    beta=atan2(yp*sin(alpha)-xp*cos(alpha),float(zp))
    return (alpha,beta)
    

def find_vector_z_rotation_transform(vect,gamma):
    """ find the rotation transform of the rotation by the euler
    angles (alpha,beta,gamma) where
    alpha and beta are determined by assuming that a vector
    (0,0,z) will be rotated to point in the direction of vect.
    See rotations.lyx for the mathematics.

    Returns a rotation matrix like find_rotation_transform
    """
    (alpha,beta)=find_vector_z_ab(vect)
    return find_rotation_transform((alpha,beta,gamma))

def find_vector_zx_rotation_abg(newZvect,newXvect):
    """
    Find the rotation transform of the rotation by the euler angles
    (alpha,beta,gamma) where the old z-axis ends up along the
    new zvect and the old x-axis end up along the new xvect. Basically
    a vector in z is now along newzvect and a vector along x is now
    along newxvect.

    see SHID/Documents/MathTricks/Coordinate Systems/rotations.lyx
    for the mathematics
    """
    
    
    xAp=float(newZvect.x())
    yAp=float(newZvect.y())
    zAp=float(newZvect.z())
    z=float(newZvect.magnitude())

    xBp=float(-newXvect.x())
    yBp=float(-newXvect.y())
    zBp=float(-newXvect.z())
    alpha=atan2(-yAp,float(xAp))
    beta=atan2(yAp*sin(alpha)-xAp*cos(alpha),float(zAp))

    if xAp == 0.0:
        gamma=atan2(yBp,xBp)       
    elif xBp == 0.0:
        gamma=atan2(yAp,-xAp)        
    else:        
        gamma=atan2((zAp/float(z))*((yAp/float(xAp))-(yBp/float(xBp))),1+(yBp*yAp/float(xBp*xAp)))
    return (alpha,beta,gamma)
    
def find_vector_zx_rotation_transform(newZvect,newXvect):

    return find_rotation_transform(find_vector_zx_rotation_abg(newZvect,newXvect))


def _test_find_vector_zx_rotation_transform():
    from random import uniform
    a2=uniform(0,2*pi)
    b2=uniform(0,2*pi)
    g2=uniform(0,2*pi)
    eulerangles=[a2,b2,g2]
    #print "secret euler angles:",[e*180/pi for e in eulerangles]
    rotmat_body_to_sensor=find_rotation_transform(eulerangles)

    vz=Vector((0,0,1),vtype="CARTESIAN")
    vx=Vector((1,0,0),vtype="CARTESIAN")

    vz_measured=vz.rotate(rotmat_body_to_sensor)
    vx_measured=vx.rotate(rotmat_body_to_sensor)

    eulerangles_found=list(find_vector_zx_rotation_abg(vz_measured,vx_measured))
    #print "found eulerangles:",[e*180/pi for e in eulerangles_found]
    
    rotmat_sensor_to_body=find_rotation_transform(eulerangles_found)
    rotmat_sensor_to_body=invert(array(rotmat_sensor_to_body))

    vz_compensated=vz_measured.rotate(rotmat_sensor_to_body)
    vx_compensated=vx_measured.rotate(rotmat_sensor_to_body)
    #print "vzmeasured:",vz_measured
    print "vzcompenasted:",vz_compensated

    #print "vxmeasured:",vx_measured
    print "vxcompenasted:",vx_compensated

        

class Triangle(BaseClass):
    """ 
    three points in space forming a triangle. The points are ordered according
    to the right hand rule so that the triangle normal points in the direction
    of a right handed thumb if your fingers curl in the direction of the point
    order.
    """
    def __init__(self,points):
        if len(points)!=3:
            raise Exception
        attributesFromDict(locals())
    def normalVect(self):
        # I am not certain that this is right but here goes...
        return self.points[0].cross(self.points[1]) 

class PointGrid(BaseClass):
    """ This is a rectangular grid of points """
    def __init__(self,xstart,xend,nx,ystart,yend,ny,zstart,zend,nz):
        import scipy
        self.points=[]
        self.nx=nx
        self.ny=ny
        self.nz=nz
        if xstart==xend:
            dx=0
            self.xvals=[xstart]
        else:
            dx=(xend-xstart)/float(nx-1)
            self.xvals=scipy.arange(xstart,xend+dx/2.0,dx)
        if ystart==yend:
            dy=0
            self.yvals=[ystart]
        else:
            dy=(yend-ystart)/float(ny-1)
            self.yvals=scipy.arange(ystart,yend+dy/2.0,dy)
        if zstart==zend:
            dz=0
            self.zvals=[zstart]
        else:
            dz=(zend-zstart)/float(nz-1)
            self.zvals=scipy.arange(zstart,zend+dz/2.0,dz)
        for x in self.xvals:
            for y in self.yvals:
                for z in self.zvals:
                    self.points.append(Vector((x,y,z)))


class Shape(BaseClass):
    pass


class Sphere(BaseClass):
    def __init__(self,posvect,radius,_ClassName_="Sphere"):
        self._ClassName_=_ClassName_
        if not isinstance(posvect,Vector):
            raise Exception
        self.Pos=posvect
        self.Radius=radius

class LineSeg(Shape):
    """ Defines a Line Segment in terms of 'points'"""

    XVector=Vector((1.0,0.0,0.0))
    YVector=Vector((0.0,1.0,0.0))
    ZVector=Vector((0.0,0.0,1.0))

    def __init__(self,startvect,endvect,_ClassName_="LineSeg"):
        self._ClassName_=_ClassName_
        if not isinstance(startvect,Vector):
            raise Exception
        if not isinstance(endvect,Vector):
            raise Exception
        self.startvect=startvect
        self.endvect=endvect
        self._center = None
        self._length = None
        self._vector = None
        self._unitvector = None
        # See SHID/Documents/MathTricks/VectorsandField.lyx for the line parameters
        self.a=self.endvect.x()-self.startvect.x()
        self.b=self.endvect.y()-self.startvect.y()
        self.c=self.endvect.z()-self.startvect.z()
        return

    def xyzStartEnd(self):
        return [self.startvect.xyz(),self.endvect.xyz()]


    def two_orthonorm(self,verbose=False):
        """
        This function produces two 
        unit vectors that are orthogonal
        to this line segment
        """
        v=self.unitvector()
        diffv=v-LineSeg.ZVector
        
        if basicallyzero(diffv.x()) and basicallyzero(diffv.y()):
            #Our vector points along z
            uv1=LineSeg.XVector
            uv2=LineSeg.YVector        
        else:
            uv1=v%LineSeg.ZVector
            uv1=uv1.unitvector()
            uv2=v%uv1
            uv2=uv2.unitvector()
        if verbose:
            
            print "UV1:",uv1
            print "UV2:",uv2
        return uv1,uv2
        
        

    def center(self):
        if self._center == None:
            self._center = self.compute_center()
        return self._center

    def length(self):
        if self._length == None:
            self._length,self._vector=self.compute_length_direction()
        return self._length

    def vector(self):
        if self._vector == None:
            self._length,self._vector=self.compute_length_direction()
        return self._vector
            

    def unitvector(self):
        if self._unitvector == None:
            self._unitvector=self.vector().unitvector()
        return self._unitvector


    def __str__(self,Mode="CARTESIAN"):
        return "Start: "+self.startvect.__str__(Mode=Mode)+" End: "+self.endvect.__str__(Mode=Mode)

    def shortestDistance(self,vect,return_PerpLineSeg=False):
        """
        As discussed in SHID/Documents/MathTricks/VectorsandFields.lyx
        this function returns the shortest distance from the point described
        by vect and the LineSeg. If the return_PerpLineSeg==True then
        the function returns a Line Segment starting at the vect and
        ending at the point closest to vect on the line.
        """
        denom=self.a*self.a+self.b*self.b+self.c*self.c
        t=(self.a*(vect.x()-self.startvect.x())+self.b*(vect.y()-self.startvect.y())+self.c*(vect.z()-self.startvect.z()))/float(denom)
        closestpoint=Vector((self.startvect.x()+self.a*t,
                             self.startvect.y()+self.b*t,
                             self.startvect.z()+self.c*t),type="CARTESIAN")

        perplineseg=LineSeg(vect,closestpoint)
        if return_PerpLineSeg==True:
            return perplineseg
        else:
            return perplineseg.length

    def negate(self):
        return Vector((-self.x,-self.y,-self.z),type="CARTESIAN")

    def translate(self,vect):
        ns=self.startvect.translate(vect)
        ne=self.endvect.translate(vect)
        return LineSeg(ns,ne)
        
    def rotate(self,transformmatrix):
        ns=self.startvect.rotate(transformmatrix)
        ne=self.endvect.rotate(transformmatrix)
        return LineSeg(ns,ne)


    def transform(self,orient):
        if hasattr(orient[0],'__len__'):
            # The first one is the rotation matrix
            return self.transformRT(orient[0],orient[1])
        else:
            # The first one is the translation vector
            return self.transformTR(orient[0],orient[1])


    def transformRT(self,rotmatrix,transvect):
        """ first rotate then translate """
        ns=self.startvect.transformRT(rotmatrix,transvect)
        ne=self.endvect.transformRT(rotmatrix,transvect)
        return LineSeg(ns,ne)
    
    def transformTR(self,transvect,rotmatrix):
        """ first translate then rotate """
        ns=self.startvect.transformTR(transvect,rotmatrix)
        ne=self.endvect.transformTR(transvect,rotmatrix)
        return LineSeg(ns,ne)

        
    def compute_center(self):
        """ Returns the line segment's center as a point"""
        xc=(self.startvect.x()+self.endvect.x())/2.0
        yc=(self.startvect.y()+self.endvect.y())/2.0
        zc=(self.startvect.z()+self.endvect.z())/2.0
        return Vector((xc,yc,zc),vtype="CARTESIAN")

    def compute_length_direction(self):
        """ Returns the line segment's length"""
        dx=self.endvect.x()-self.startvect.x()
        dy=self.endvect.y()-self.startvect.y()
        dz=self.endvect.z()-self.startvect.z()
        l=sqrt(dx*dx+dy*dy+dz*dz)
        v=Vector((dx,dy,dz),vtype="CARTESIAN")
        return l,v

    def MutualParams(self,otherlineseg):
        """ Assuming that the current line seg goes from A to B and the other
        line seg goes from a to b then this function spits out the following
        lengths:
        l=AB
        m=ab
        R1=Bb
        R2=Ba
        R3=Aa
        R4=Ab
        as described on page 55 of Grover's inductance computations
        """
        l=self.length()
        m=otherlineseg.length()
        A=self.startvect
        B=self.endvect
        a=otherlineseg.startvect
        b=otherlineseg.endvect
        R1=LineSeg(B,b).length()
        R2=LineSeg(B,a).length()
        R3=LineSeg(A,a).length()
        R4=LineSeg(A,b).length()
        return (l,m,R1,R2,R3,R4)
        
   
    def rotatePoint(self,point,angleRadians):
        """
        Given a point described with a vector,
        this function returns a new vector which
        results from rotating that point about this
        line.

        see http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/
        """
        x,y,z=point.xyz()
        a,b,c=self.startvect.xyz()
        T=angleRadians
        u,v,w=self.vector().unitvector().xyz()
        X=(a*(v**2+w**2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(T))+x*cos(T)+(-c*v+b*w-w*y+v*z)*sin(T)
        Y=(b*(u**2+w**2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(T))+y*cos(T)+(c*u-a*w+w*x-u*z)*sin(T)
        Z=(c*(u**2+v**2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(T))+z*cos(T)+(-b*u+a*v-v*x+u*y)*sin(T)
        return Vector((X,Y,Z),vtype="CARTESIAN")
        

    def reverse(self):
        """
        return a reversed copy of ourselves
        """
        return LineSeg(self.endvect,self.startvect)

class HelixSeg(LineSeg):
    """
    The helix segment is defined by two points on the surface of
    a cylinder. The helix segment does not trace out more than
    2pi of angle around the cylinder and always travels counter-clockwise
    as seen from looking down from the endvect of the cylinder (as opposed
    to the startvect of the cylinder)

    """
    def __init__(self,startvect,endvect,cylind,Mode="NORMAL"):
        if Mode=="NORMAL":
            from gsl import fcmp
            # We need to ensure that startvect and endvect
            # both point to positions on the cylinder.
            Ra=cylind.shortestDistance(startvect,return_PerpLineSeg=True)
            Rb=cylind.shortestDistance(endvect,return_PerpLineSeg=True)
            if not (fcmp(Ra.length(),cylind.radius,1e-13)==0 and fcmp(Rb.length(),cylind.radius,1e-13)==0):
                print "Endpoints are not on the cylinder"
                raise Exception
            self.cylind=cylind
            LineSeg.__init__(self,startvect,endvect)
            # The cylind_startvect and cylind_endvect are vectors described in a
            # coordinate system where x points along Ra, z points along the cylinder
            # axis and the origin is at cylind.startvect
            zcylindstart_vect=Ra.endvect.subtract(self.cylind.startvect)# points from the z position on the cylinder of helix start to cylindstart
            zcylindend_vect=Rb.endvect.subtract(self.cylind.startvect)  # points from the z position on the cylinder of helix end to cylindstart
            if zcylindstart_vect.unitvector() == None:
                zcylindstart=0
            else:
                zcylindstart=zcylindstart_vect.magnitude*((zcylindstart_vect.unitvector()).dot(self.cylind.vector.unitvector()))
            if zcylindend_vect.unitvector() == None:
                zcylindend=0
            else:
                zcylindend=zcylindend_vect.magnitude()*((zcylindend_vect.unitvector()).dot((self.cylind.vector()).unitvector()))
            crossProd=(Ra.vector()).cross(Rb.vector()).dot((self.cylind.vector()).unitvector())
            dotProd=(Ra.vector()).dot(Rb.vector())
            phi=atan2(crossProd,dotProd)
            #if phi<0.0:
            #    phi=phi+2*pi
            if phi==0.0 or phi==fabs(2*pi):
                print "phi:",phi
                raise Exception
            self.cylind_startvect=Vector((Ra.length,0,zcylindstart),vtype="CYLINDRICAL")
            self.cylind_endvect=Vector((Rb.length,phi,zcylindend),vtype="CYLINDRICAL")

            self.M=(self.cylind_endvect.z-self.cylind_startvect.z)/float(self.cylind_endvect.phi-self.cylind_startvect.phi)
            self.forward_transform=find_vector_zx_rotation_transform(self.cylind.vector,Ra.vector.negate())
            import numpy.linalg
            self.inverse_transform=numpy.linalg.inv(self.forward_transform)
            print "Phi END:",phi
            
            
        elif Mode=="CYLINDER_COORDS":
            # startvect and endvect are in coordinates referenced from cylind
            #
            print "I have not dealt with this yet"
            raise Exception
        else:
            #Unknown Mode
            print "Unknown Mode"
            raise Exception

    def toLineSegs(self,min_dphi):
        """
        Break the helix segment into a number 
        of line segments.
        """
        result=None
        phis=startvect.phi(format="POS")
        phie=endvect.phi(format="POS")
        phidiff=phie-phis
        if abs(phidiff)>=180:
            #
            print  "I think we wrapped around 360"
            print  "I have not considered this yet"
            raise Exception
        nsteps=abs(phidiff)/float(min_dphi)
        
        if nsteps<=1:
            result=[LineSeg(self.startvect,self.endvect)]
        else:
            nsteps=int(nsteps)
            zdiff=endvect.z()-starvect.z()
            dphi=phidiff/float(nsteps)
            dz=zdiff/float(nsteps)

            result=[]
            svect=self.startvect
            for i in range(nsteps):
                evect=Vector((svect.rcylind(),
                              svect.phi()+dphi,
                              svect.z()+dz),
                             vtype="CYLINDRICAL")
                tmp=LineSeg(svect,evect)
                svect=evect
                result.append()

        return result
            
            


    def transform_from_helixcoords(self,tvect):
        """
        transforms a vector from the helix-segment coordinate system
        to the coordinate system originally used to define the helix-segment
        and its cylinder.        
        """
        rotmatrix=self.inverse_transform
        return tvect.transformRT(rotmatrix,self.cylind.startvect)
        #print "Not Yet Ready"
        #raise Exception

    def transform_to_helixcoords(self,tvect):
        """
        transforms a vector to the helix-segment coordinate system
        from the coordinate system originally used to define the helix-segment
        and its cylinder.
        """
        rotmatrix=self.forward_transform
        return tvect.transformTR(self.cylind.startvect.negate(),rotmatrix)
        #print "Not Yet Ready"
        #raise Exception

    def __str__(self,AbsRel="RELATIVE",Mode="CYLINDRICAL"):
        if AbsRel=="RELATIVE":
            return "Relative to Cylinder -- Start: "+self.cylind_startvect.__str__(Mode=Mode)+" End: "+self.cylind_endvect.__str__(Mode=Mode)
        elif AbsRel=="ABSOLUTE":
            return LineSeg.__str__(self)
        else:
            print "Unknown AbsRel"
            raise Exception

class Arc(Shape):
    """

    """
    def __init__(self,
                 center,
                 start_angle,
                 end_angle,
                 radius,
                 axis=Vector((0,0,1)),
                 _ClassName_="Arc",
                 ):
        attributesFromDict(locals())
        # our job is now to figure out the startpos and endpos
        lseg=LineSeg(center,center+axis)
        r0,rperp=lseg.two_orthonorm()
        r0=r0.scale(radius)
        if (axis-Vector((0,0,1))).length() >= 1e-9:
            raise Exception
        
        rotA=find_yaw_transform(start_angle)
        #rotA=find_rotation_about_axis(axis,start_angle)
        print "rotA",rotA
        self.startvect=center+r0.rotate(rotA)

        rotB=find_yaw_transform(end_angle)
        #rotB=find_rotation_about_axis(axis,end_angle)
        self.endvect=center+r0.rotate(rotB)
        

    def reverse(self):
        """
        return a reversed copy of ourselves
        """
        return Arc(self.center,
                   self.end_angle,
                   self.start_angle,
                   self.radius,
                   axis=self.axis)
        
class RectPrisim(Shape):
    """
    width is perpendicular to both lvect and hvect
    """
    def __init__(self,centervect,lvect,hvect,wvect):
        #wvect=lvect.cross(hvect).unitvector()
        #wvect=wvect.scale(width)
        length=lvect.length()
        height=hvect.length()
        width=wvect.length()
        attributesFromDict(locals())

    def NodePoints(self):
        """
        give a list of the eight nodal points
        """
        Point1=self.centervect+self.lvect.scale( 0.5)+self.hvect.scale( 0.5)+self.wvect.scale( 0.5)
        Point2=self.centervect+self.lvect.scale(-0.5)+self.hvect.scale( 0.5)+self.wvect.scale( 0.5)
        Point3=self.centervect+self.lvect.scale(-0.5)+self.hvect.scale(-0.5)+self.wvect.scale( 0.5)
        Point4=self.centervect+self.lvect.scale( 0.5)+self.hvect.scale(-0.5)+self.wvect.scale( 0.5)
        Point5=self.centervect+self.lvect.scale( 0.5)+self.hvect.scale( 0.5)+self.wvect.scale(-0.5)
        Point6=self.centervect+self.lvect.scale(-0.5)+self.hvect.scale( 0.5)+self.wvect.scale(-0.5)
        Point7=self.centervect+self.lvect.scale(-0.5)+self.hvect.scale(-0.5)+self.wvect.scale(-0.5)
        Point8=self.centervect+self.lvect.scale( 0.5)+self.hvect.scale(-0.5)+self.wvect.scale(-0.5)
        return [Point1,Point2,Point3,Point4,Point5,Point6,Point7,Point8]
        
class Cylinder(LineSeg):
    """ Defines a Cylinder Segment in terms of 'vectors'"""
    def __init__(self,startpoint,endpoint,radius,_ClassName_="Cylinder"):
        self.radius=radius
        LineSeg.__init__(self,startpoint,endpoint,_ClassName_=_ClassName_)

class CylindericalShell(Cylinder):
    """ Defines a Cylinderical Shell in terms of 'points'"""
    def __init__(self,startpoint,endpoint,radius,thickness,_ClassName_="CylindericalShell"):
        self.thickness=thickness
        Cylinder.__init__(self,startpoint,endpoint,radius,_ClassName_=_ClassName_)

class Tube(Shape):
    """
    A list of vectors describing the path of the tube and some way
    to make the tube have thickness along the path.
    """
    pass
class CylindricalTube(Tube):
    def __init__(self,pathList,radii,
                 isALoop=False):
        """
        A list of vector points that the tube
        traverses and the radii that it has at
        those points. The orientation of the 
        circles of the tube must be such that
        they are the average of the two orientations
        of the lineseg on either path leading into
        the section. If this is not a closed loop
        then the start and end circles are just
        oriented with the one segment with which
        they are associated.

        If this is a loop, we close the tube back
        at the starting point.
        """
        # We need to figure out the line segments
        # and the direction that the circle points.
        attributesFromDict(locals())
        self.segs=[]
        for i in range(len(pathList)-1):
            A=pathList[i]
            B=pathList[i+1]
            self.segs.append(LineSeg(A,B))
        # now the directions at a point
        self.pointDirections=[self.segs[0].unitvector()]
        for i in range(len(self.segs)-1): # skip last segment
            A=self.segs[i].unitvector()
            B=self.segs[i+1].uintvector()
            avg=(A+B).scale(0.5)
            self.pointDirections.append(avg)
        self.pointDirections.append(self.segs[-1].unitvector())
    def GetTriangles(self,numPhi=6):
        """
        We can get a representation in triangles of this tube.
        """
        triangles=[]
        
        return triangles
        
    

class Cone(LineSeg):
    """ Defines a Cone in terms of 'vectors'"""
    def __init__(self,basecentre,heightvect,radius,_ClassName_="Cone"):
        self.radius=radius
        LineSeg.__init__(self,basecentre,basecentre+heightvect,_ClassName_=_ClassName_)

class Arrow(LineSeg):
    """ Defines a cylinder capped with a cone  """
    def __init__(self,startpoint,endpoint,cylradius,
                 coneradiusFactor=2.0,
                 coneheightFactor=2.0, #This is what you do 
                                       #if you don't want to
                                       #specify coneheight
                 _ClassName_="Arrow"):
        
        self.cylradius=cylradius
        self.coneradius=cylradius*coneradiusFactor
        self.coneheight=self.coneradius*coneheightFactor
        
        self.height=(endpoint-startpoint).length()
        
        LineSeg.__init__(self,startpoint,endpoint,_ClassName_=_ClassName_)
        self.breakpos=startpoint+self.unitvector().scale(self.height-self.coneheight)  #Where the cylinder ends and the cone starts
        #print "SP:",startpoint.xyz()
        #print "EP:",endpoint.xyz()
        #print "BP:",self.breakpos.xyz()

class Path(Shape):
    """
    Really its just a list of
    positions from first to last
    """
    def __init__(self,vects):
        if(len(vects)<2):
            raise Exception
        attributesFromDict(locals())

class Hexahedron(Shape):
    """
    page 274 of vtk book

    You may want to look at KVTK if you are making
    a lot of these. There it shows how you defines
    points first and then define hexahedrons with
    the point numbers
    """
    def __init__(self,P0,P1,P2,P3,P4,P5,P6,P7):
        #attributesFromDict(locals())
        self._P0=P0
        self._P1=P1
        self._P2=P2
        self._P3=P3
        self._P4=P4
        self._P5=P5
        self._P6=P6
        self._P7=P7

    def P0(self):
        return self._P0

    def P1(self):
        return self._P1

    def P2(self):
        return self._P2

    def P3(self):
        return self._P3

    def P4(self):
        return self._P4

    def P5(self):
        return self._P5

    def P6(self):
        return self._P6

    def P7(self):
        return self._P7
       

class Quadrilateral(Shape):
    
    def __init__(self,centerVector,areaVector,lengthVector):
        """
        centerVector is a vector describing where the center is

        lengthVect is a vector describing the direction of the long axis
        and its length

        areaVector is a vector describing the 'direction' of the
        area
        

        N1-N4 count counter-clockwise about the areaVector
        """
        attributesFromDict(locals())
        self._N1=None
        self._N2=None
        self._N3=None
        self._N4=None
        
        self.widthVector=(self.areaVector%(self.lengthVector.unitvector())).scale(1.0/float(self.lengthVector.length()))
        
    
    def N1(self):
        if self._N1==None:
            self._N1=self.centerVector+self.lengthVector.scale(-0.5)+self.widthVector.scale(0.5)
        return self._N1

    def N2(self):
        if self._N2==None:
            self._N2=self.centerVector+self.lengthVector.scale(0.5)+self.widthVector.scale(0.5)
        return self._N2

    def N3(self):
        if self._N3==None:
            self._N3=self.centerVector+self.lengthVector.scale(0.5)+self.widthVector.scale(-0.5)
        return self._N3

    def N4(self):
        if self._N4==None:
            self._N4=self.centerVector+self.lengthVector.scale(-0.5)+self.widthVector.scale(-0.5)
        return self._N4

    def area(self):
        return self.areaVector.length()

def lineseg_tests():
    l=LineSeg(Vector((0,0,0),vtype="CARTESIAN"),Vector((5,5,0),vtype="CARTESIAN"))
    P=Vector((1.5,3.5,0),vtype="CARTESIAN")
    s=l.shortestDistance(P,return_PerpLineSeg=True)
    print s," Length: ",s.length


def helixseg_tests():
    cylindA=Cylinder(Vector((0,0,0),vtype="CARTESIAN"),Vector((0,0,5),vtype="CARTESIAN"),5)
    SegA=HelixSeg(Vector((5,0,0),vtype="CARTESIAN"),Vector((0,5,2),vtype="CARTESIAN"),cylindA)
    print SegA

class Rect2D(Quadrilateral):
    def __init__(self,centerVector,lengthX,widthY):
        if centerVector.z()!=0.0:
            raise Exception
        area=lengthX*widthY
        areaVector=Vector((0,0,area),vtype="CARTESIAN")
        lengthVector=Vector((lengthX,0,0),vtype="CARTESIAN")
        Quadrilateral.__init__(self,centerVector,areaVector,lengthVector)

if __name__ == "__main__":
    #helixseg_tests()
    #lineseg_tests()
    #P=PointGrid(1,5,10,0,0,1,0,0,1)
    #Arc(Vector((8.4672,66.3125,0)),177.44*pi/180.0,358.57*pi/180.0,1.008200)
    _test_find_vector_zx_rotation_transform()
