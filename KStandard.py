#!/usr/bin/env python
#     Title:      KStandard
#     Purpose: a place to store some generally useful code.
#              



#*********** TODO *************


#********* Event Log **********
# Added identity matrix
# November 21,2006 (Karl Edler)
#-----------------------------


from __future__ import division
import cPickle
from math import fabs,pi,modf,floor,ceil,sqrt,atan2,sin,cos,acos
from pprint import pprint
from types import *


def SkimDefaults(self,defaults,kwargs):
    for d in defaults:
        if kwargs.has_key(d):
            exec("self."+d+"="+repr(kwargs[d]))
            del kwargs[d]
        else:
            exec("self."+d+"="+repr(defaults[d]))            

def DictBuilder(keys,default):
    result={}
    for k in keys:
        result[k]=default
    return result

def KBisect(lister,val,compare=cmp):
    """
    return index if found.
    otherwise return -indexafter which may be longer than the list if
    val would fit after all the values in the list.
    if val would fit before all the values in the list the result is
    -(nmemb+1)
    """
    nmemb=len(lister)
    high=nmemb-1
    low=0
    while(low<=high):
        mid=low+((high-low)/2);
        #print "mid:",mid
        comp=compare(lister[mid],val)
        if comp<0:
            # lister[mid] < val
            low=mid+1
        elif comp >0:
            # lister[mid] > val
            high=mid-1
        else:
            result=mid # val found
            return result
    if(low==0):
        # all the values in lister are
        # greater than val
        result=-(nmemb+1)
        return result
    else:
        result=-low # val not found
                    # return the index after 
                    # except negative
        return result 


def DataSort(xvals,data):
    """
    we have data associated with
    the locations in xvals but
    unfortunately the xvals are
    not in order.

    This function sorts the xvals list
    and sorts the data list at the same
    time

    """
    if len(xvals)!=len(data):
        print "Data must be the same length"
        raise Exception
    
    class SElem:
        def __init__(self,x,v):
            self.x=x
            self.v=v
            
        def __cmp__(self,other):
            return cmp(self.x,other.x)
    L=[ SElem(x,v) for x,v in zip(xvals,data) ]
    L.sort()
    
    rxvals=[l.x for l in L]
    rdata=[l.v for l in L]
    return rxvals,rdata


## def allcombos(list):
##     """returns a list [[a,b],[a,b]...] containing all the combinations in list """
##     result=[]
##     for i in range(len(list)):
##         elem=list[i]
##         smaller=list[i+1:]
##         for j in smaller:
##             result.append([elem,j])
##     return result

def listadd(a,b):
    if len(a)!=len(b):
        raise Exception
    result=[va+vb for va,vb in zip(a,b)]
    return result


def listintegrate(list,dx):
    """
    Use simpson's rule to work
    out the integral of a numerical 
    list [y0,y1,y2,...] seperated by dx.

    The resulting integral list is:
    i0=0
    i1=dx*y0+0.5(y1-y0)
    
    in=dx*0.5*(yn+yn-1)
    """
    sum=0
    result=[sum]
    vprev=list[0]
    for v in list[1:]:
        contribution=dx*0.5*(v+vprev)
        sum=sum+contribution
        result.append(sum)
        vprev=v
    return result
        
def RunningAverage(lister):
    total=0
    result=[]
    for i,v in enumerate(lister):
        num=i+1
        total=total+v
        result.append(total/float(num))
    return result
        



def numberpadder(numchars,number):
    number=number.__str__()
    diff=numchars-len(number)
    if diff<0:
        raise Exception
    if diff>0:
        prepend="0"
        for i in range(diff-1):
            prepend=prepend+"0"
        number=prepend+number
    return number
        




def identitymatrix(n):
    """ Produces a square identity matrix of rank n"""
    result=[]
    for i in range(n):
        temp=[]
        for j in range(n):
            if j==i:
                temp.append(1)
            else:
                temp.append(0)
        result.append(temp)
    return result

#def xyz_from_rtp_vect(Vrtp,rtp):
#    """Converts a vector of format r,theta,phi to format x,y,z"""
#    r,t,p=rtp
#    x_hat,y_hat,z_hat=spherical_unitvectors(t,p)
#    Vx=dot(Vrtp,x_hat)
#    Vy=dot(Vrtp,y_hat)
#    Vz=dot(Vrtp,z_hat)
#    return [Vx,Vy,Vz]
    

def xyz_to_rtp_point(xyz):
    """Converts the point coordinates xyz=[x,y,z] to point coordinates
    rtp=[r,theta,phi] in spherical coordinates"""
    x=xyz[0]
    y=xyz[1]
    z=xyz[2]
    rtp=[sqrt(x*x+y*y+z*z),atan2(sqrt(x*x+y*y),z),atan2(y,x)]
    return rtp

def spherical_unitvectors(theta,phi):
    """Given the spherical coordinates theta and phi from the back of
    Griffiths this function returns the components of the cartesian unit
    vectors in terms of the spherical unitvectors.

    The result is:[[r_hat,theta_hat,phi_hat],[r_hat,theta_hat,phi_hat],[r_hat,theta_hat,phi_hat]]

    Which corresponds to [x_hat,y_hat,z_hat]
    """
    x_hat=[sin(theta)*cos(phi),cos(theta)*cos(phi),-sin(phi)]
    y_hat=[sin(theta)*sin(phi),cos(theta)*sin(phi),cos(phi)]
    z_hat=[cos(theta),-sin(theta),0.0]
    return [x_hat,y_hat,z_hat]

def allcombos(list):
    """the result will be a list of pairs of elements that were in the 'list'
    Every possible pair will be in the result. An element will not be paired
    with itself."""
    result=[]
    for i,elem in enumerate(list):
        for o in list[i+1:]:
            result.append([elem,o])
    return result

def kround(number):
    frac,whole=modf(number)
    if frac < 0.50:
        return int(floor(number))
    else:
        return int(ceil(number))

def is_odd(number):
    n=number%2
    if n == 1:
        return True
    else:
        return False

CCW=0
CW=1

def arcdirection(seg):
    """Returns the direction of the rotation as either CW or CCW"""
    xs=seg[1][0]
    ys=seg[1][1]
    xe=seg[2][0]
    ye=seg[2][1]
    diff=xs*ye-ys*xe
    if diff>0: return CCW
    else: return CW

def arcdirectionold(seg):
    """Returns the direction of the rotation as either CW or CCW"""
    xs=seg[1][0]
    ys=seg[1][1]
    xe=seg[2][0]
    ye=seg[2][1]
    ts=atan2(ys,xs)
    te=atan2(ye,xe)

    if ts<0: ts=ts+2*pi
    if te<0: te=te+2*pi


    if abs(ts-te)==pi: raise Exception
    diff=te-ts
    if diff>0:
        if pi>diff:
            return CCW
        else: return CW
    else:
        diff=-diff
        if pi>diff:
            return CW
        else:
            return CCW

    raise Exception

def circlecloser(theta,theta1,theta2):
    """Given all three values are on the range 2pi>theta>0 or on the range
    pi>theta>-pi and a theta to test and two other thetas this function
    will return:
    True if theta is closer to theta1 and
    False if theta is closer to theta2"""
    if theta<0: theta=theta+2*pi
    if theta1<0: theta1=theta1+2*pi
    if theta2<0: theta2=theta2+2*pi

    diff1=abs(theta-theta1)
    diff2=abs(theta-theta2)

    if diff1<diff2: return True
    else: return False
    
def circlecompare(theta,thetamax,thetamin):
    """Given all three values are on the range 2pi>theta>0 or on the range
    pi>theta>-pi this function will return true
    if rotating counter clockwise from thetamin to thetamax encounters
    theta even if the counterclockwise rotation passes through zero

    returns true if theta is encountered
    otherwise returns false
    """
    if theta<0: theta=theta+2*pi
    if thetamax<0: thetamax=thetamax+2*pi
    if thetamin<0: thetamin=thetamin+2*pi

        
    if thetamin>thetamax: thetamax=thetamax+2*pi
    if thetamin>theta: theta=theta+2*pi
    if thetamax>theta>thetamin:
        return True
    else:
        return False


def ListSum(numberlist):
    """Given a list of numbers this function returns the sum"""
    result = 0
    for num in numberlist: result = result + num
    return result

def basicallyzero(val,epsilon=1e-13):
    if(val<epsilon and val>-epsilon):
        return True
    else:
        return False

def fgeq(x,y,epsilon=1e-13):
    from gsl import fcmp
    "Returns true if x is greater than or equal to y"
    val=fcmp(x,y,epsilon)
    if val>=0:
        return True
    else:
        return False


def fleq(x,y,epsilon=1e-13):
    from gsl import fcmp
    "Returns true if x is less than or equal to y"
    val=fcmp(x,y,epsilon)
    if val <=0:
        return True
    else:
        return False


def angular_pos_normalizer(radians):
    """Given a position in radians this function converts it to a value between
    0 and 2pi"""
    by2pi=radians/(2*pi)
    frac,int=modf(by2pi)
    radians=2*pi*frac    
    if radians < 0:
        radians = radians + 2*pi
    return radians

def is_even(number):
    return not (is_odd(number))

def dictFromNamedObjects(List):
    result={}
    for o in List:
        result[o.name]=o
    return result

def attributesFromDict(d):
    """
    Automatically instantiate member variables
    from constructor arguments. It works like this:

    class foo:
       def __init__(self,A,B,C):
          " self.A,self.B,self.C are all set "           
          attributesFromDict(locals())

    """
    self=d.pop('self')
    for n,v in d.iteritems():
        setattr(self,n,v)

class BaseClass(object):
    """This Class should be the Great Grand-Daddy of all classes that
    I define. It will have store whatever neat tricks I think that all
    classes should have"""
    _ClassName_="BaseClass"#Introspective code should be able to tell the class name
    _Name=None             # If this is a named object we should be able to tell.
    _cPickleNotXML_=[]     # a list of sequences (not normal attributes) that should be saved and loaded with cPickle not XML
   
   

    def __init__(self):
        """It is important that this function be callable without arguments to
        make an empty instance. This can be done by setting default arguments"""
        pass

##     def __initempty__(self):
##         """Initialize everything except the member variables.
##         Setting the defaults is okay"""
##         BaseClass.__init__(self)


       

    def XML_Rebuild_Object_Top(self,XML_Pickle_Load_Object,execspace):
        """This top level ensures that we are loading the correct kind of class"""
        if not hasattr(XML_Pickle_Load_Object,"_ClassName_"):
            print "Every one of my XML files should correspond to an object"
            raise Exception
        if self._ClassName_ != XML_Pickle_Load_Object._ClassName_:
            print "The top level class should fit!"
            raise Exception
        return self.XML_Rebuild_Object(XML_Pickle_Load_Object,execspace)
                       

    def XML_Rebuild_Object(self,XML_Pickle_Load_Object,execspace):
        """Given an object of the form constructed by XML_Pickle_load this will
        return an object that has the correct class types"""
        # Construct a Class for this object and every sub-object that
        # has the attribute _ClassName_
        # place the data members in the new object
        xml_obj=XML_Pickle_Load_Object
        if not hasattr(xml_obj,"_ClassName_"):
            # If it doesn't have a class name and is not a list
            # I can't do much with it.
            if not isSequence(xml_obj):
                return xml_obj
            else:
                return [self.XML_Rebuild_Object(j,execspace) for j in xml_obj]
        # Make a new object that will be filled in and passed back
        klass=xml_obj._ClassName_
        #print klass
        exec "result="+klass+"()" in execspace,locals()
        # Develop a list of data attributes of xml_obj
        dataobjects = [(o,getattr(xml_obj,o)) for o in dir(xml_obj) if not callable(getattr(xml_obj,o))]
        #pprint(dataobjects)
        for dname,d_xml in dataobjects:
            #print dname
            d=self.XML_Rebuild_Object(d_xml,execspace)
            exec "result."+dname+"=d"
        return result

    def dump(self,filename):
        """This method dumps the contents of the class to a file"""
        fh=open(filename,"w+")
        cPickle.dump(self,fh)
        fh.close()

    def dumps(self):
        return cPickle.dumps(self)

    def load(self,filename):
        """This method returns the object contained in a file"""
        fh=open(filename,"r")
        newobject=cPickle.load(fh)
        fh.close()
        #print self
        #print newobject
        #self=newobject
        return newobject

    def loads(self,string):
        return cPickle.loads(string)

    def XML_load(self,filename,execspace):
        """
        newobj=classtype().XML_load(filename,globals())
        """
        xml_obj=self.XML_Pickle_load(filename)
        return self.XML_Rebuild_Object_Top(xml_obj,execspace)

    def XML_Objectify_load(self,filename):
        stream=file(filename,"r")
        xml_obj=XML_Objectify(filename)
        stream.close()
        return xml_obj.make_instance()
        


    def XML_dump(self,filename):
        if self._ClassName_ == None:
            print "Will Be unable to unpickle XML properly"
            raise Exception
        stream=file(filename,"w+")
        XML_Pickler(self).dump(stream)
        stream.close()

    def XML_Pickle_load(self,filename):
        stream=file(filename,"r")
        result=XML_Pickler().load(stream)
        stream.close()
        return result




    def __str__(self):
        """We need some way to print ourselves out properly"""
        if not self._Name == None:
            return self._Name
        else:
            return "UnNamed"#str(self)
    




def getdatestring():
    import commands
    output=commands.getoutput("date")
    return output

def sign(val):
    if val > 0: return 1
    elif val < 0: return -1
    else: return 0

def listaverage(list):
    num=len(list)
    result=sum(list)/float(num)
    return result


def sign_majority_average(points):
    """Given a list of values this function finds the sign of the majority. Then it finds the
    average of all the values with this sign and returns that average."""
    negatives=[]
    positives=[]
    if is_even(len(points)):
        print "There should be an odd number of points"
        raise Exception
    for p in points:
        if p > 0:
            positives.append(p)
        else:
            negatives.append(p)

    nump=len(positives)
    numn=len(negatives)
    if nump > numn:
        result=listaverage(positives)
    else:
        result=listaverage(negatives)
    return result

def absmax_index(list):
    """Find the absolute value maximum in the list"""
    maxer=0
    index=0
    for i,e in enumerate(list):
        ne=fabs(e)
        if ne>maxer:
            maxer=ne
            index=i
    return [index,maxer]

def absmax(list):
    return absmax_index(list)[1]

def matrix_abs_max(matrix):
    import numpy
    result=None
    for v in matrix:
        
        #print "Entered Matrix Max"
        if type(v) is numpy.ndarray or type(v) is ListType:
            temp=matrix_abs_max(v)
            if result==None:
                result=temp
            elif temp>result:
                result=temp
        else:
            if result==None:
                result=abs(v)
            elif abs(v)>result:
                result=abs(v)
    return result                  

def matrix_abs_min(matrix):
    import numpy
    result=None
    for v in matrix:
        if type(v) is numpy.ndarray or type(v) is ListType:
            temp=matrix_abs_min(v)
            if result==None:
                result=temp
            elif temp<result:
                result=temp
        else:
            if result==None:
                result=abs(v)
            elif abs(v)<result:
                result=abs(v)
    return result

def matrix_max(matrix):
    """ Find the maximum even if we have lists of lists """
    import numpy
    result=None
    for v in matrix:
        #print "Entered Matrix Max"
        if type(v) is numpy.ndarray or type(v) is ListType:
            temp=matrix_max(v)
            if result==None:
                result=temp
            elif temp>result:
                result=temp
        else:
            if result==None:
                result=v
            elif v>result:
                result=v
    return result                  

def matrix_min(matrix):
    """ Find the manimum even if we have lists of lists """
    import numpy
    result=None
    for v in matrix:
        if type(v) is numpy.ndarray or type(v) is ListType:
            temp=matrix_min(v)
            if result==None:
                result=temp
            elif temp<result:
                result=temp
        else:
            if result==None:
                result=v
            elif v<result:
                result=v
    return result
    

def inclusive_range(start,end):
    """This function produces an integer range of values from start to end
    inclusive"""
    result=range(start,end+1)
    return result


def info(object, spacing=10, collapse=1):
    """Print methods and doc strings.

    Takes module, class, list, dictionary, or string.
    """

    methodList = [method for method in dir(object) if callable(getattr(object,method))]
    processFunc = collapse and (lambda s: " ".join(s.split())) or (lambda s: s)

    print "\n".join(["%s %s" %
                     (method.ljust(spacing),
                     processFunc(str(getattr(object,method).__doc__)))
                     for method in methodList])

def has_method(object,attr):
    methodList = [method for method in dir(object) if callable(getattr(object,method))]
    if methodList.__contains__(attr):
        return True
    else:
        return False
    

def datainfo(object,suppressHidden=True):
    """returns a dictionary of data values"""
    dataobjects={}
    for o in dir(object):
        if not callable(getattr(object,o)):
            if suppressHidden and o[:2]=="__" and o[-2:]=="__":
                pass
            else:
                dataobjects[o]=getattr(object,o)
    return dataobjects
    


def wait(str=None, prompt='Press return to continue...\n'):
    if str is not None:
        print str
    raw_input(prompt)
    
def buildConnectionString(params):
    """Build a connection string from a dictionary of paramters.


    Returns string."""
    if not isDictionary(params): raise TypeError
    return ":".join(["%s=%s" % (k,v) for k, v in params.items()])

def listmult(list,factor):
    return [factor*q for q in list]

def isDictionary(param):
    t={}
    if type(param)==type(t): return True
    else: return False
if __name__ == "__main__":
    myParams = {"server":"mpilgrim",\
                "database":"master",\
                "uid":"sa",\
                "pwd":"secret"}
    print buildConnectionString(myParams)

def isTuple(param):
    t=()
    if type(param)==type(t): return True
    else: return False

def isList(param):
    t=[]
    if type(param)==type(t): return True
    else: return False

def isSequence(param):
    if isTuple(param) or isList(param): return True
    else: return False

def flatten(sequence):

        def hms(fpd):
                if fpd < 60:
                    return fpd
                elif fpd < 60**2:
                    return "%s:%s" % (int(fpd/60), fpd-int(fpd/60)*60)
                else:
		    h = int(fpd/60**2)
		    fpd -= h*60**2
		    m = int(fpd/60)
		    fpd -= m*60
		    s = fpd
                    return "%s:%s:%s" % (h, m, s)
                
	def rflat(seq2):
		seq = []
		for entry in seq2:
			if seqin([entry]):
        			seq.extend([i for i in entry])
			else:
				seq.append(entry)
		return seq

	def seqin(sequence):
		for i in sequence:
			if ('__contains__' in dir(i) and    ## all sequences have '__contains__' in their dir()
                           type(i) != str and type(i) != dict): ## parentheses present to aid commenting mid-condition
				return True
		return False

        import time
        btime = time.time()
        d1time = btime
	seq = [sequence][:]    ## in case parameter isn't already a sequence
	print "Thinking",
	while seqin(seq):
                d2time = time.time()
                if d2time-d1time >= 5:
                    print ".",
                    d1time = d2time
		seq = rflat(seq)
	atime = time.time()
	print
	print "Sequence flattened in " + str(hms(atime-btime))
	return seq


def findclassname(q):
    return q.__class__.__name__

def matrix_max_min_tester():
    q=[1000,2,3,4,[10,2,0],0,2,[100,-100]]
    print matrix_max(q)
    print matrix_min(q)

def LineSegTester():
    AB=LineSeg(Point((1,2,3),type="CARTESIAN"),Point((2,3,1),type="CARTESIAN"))
    ab=LineSeg(Point((3,2,1),type="CARTESIAN"),Point((1,3,2),type="CARTESIAN"))
    print AB.MutualParams(ab)

if __name__ == "__main__":
    #print numberpadder(5,30.0)
    print info.__doc__ 
    print listintegrate([0,1,2,3,4,5],1)
    #LineSegTester()
    #matrix_max_min_tester()
    #p=point((0,0,0),type="CARTESIAN")
    #p.translate(point((0,0,-5.0),type="CARTESIAN"))
    #print p.x," ",p.y," ",p.z

    
