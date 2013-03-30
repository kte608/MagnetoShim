#!/usr/bin/env python
import pylab
from matplotlib import rc
import matplotlib.pyplot as plt
#rc('text',usetex=True)
import cPickle
from KStandard import is_odd
from scipy import array,zeros,arange
from pprint import pprint

def HarmonicBarChart(mchoice,harmonicResponse,showPrime=False,titleappend="",mchoiceMinusOne=True,ymax=None,ymin=None):
    datadict={}
    def keygen(n,m):
        return n.__str__()+"_"+m.__str__()
    nmax=None
    nmin=None

    for d in harmonicResponse:
        n=int(d[0])
        m=int(d[1])
        if nmin==None:
            nmin=n
        if nmax==None:
            nmax=n
        if n<nmin:
            nmin=m
        if n>nmax:
            nmax=n
        datadict[keygen(n,m)]=d[2:]
    barwidth=1.0
    if showPrime:
        q=1
        titlestring="$T^{\prime}_{n,"+mchoice.__str__()+"}$"
    else:
        q=0
        titlestring="$T_{n,"+mchoice.__str__()+"}$"
    titlestring=titlestring+titleappend#"\n  Harmonic Response: "+ResponseFile
    if mchoiceMinusOne:
        nlist=range(mchoice-1,nmax)
    else:
        nlist=range(mchoice,nmax)
    responses=[abs(datadict[keygen(n,mchoice)][q]) for n in nlist]
    #responses=responses[:7]
    #print "nlist:"
    #pprint(nlist)
    #pprint(responses)
    num=len(responses)
    indicies=arange(num)
    #print "num:",num
    #pprint(indicies)
    ax=plt.subplot(111)
    
    
    plt.title(titlestring)
    plt.ylabel("response")
    plt.xlabel("n")
    
    
    plt.xticks(indicies+barwidth*0.5,[i.__str__() for i in nlist])
        
    ax.set_yscale("log")
    plt.bar(indicies,responses,barwidth)
    #plt.axis([min(indicies),max(indicies),min(responses)*0.1,1])
    if ymax==None:
        ymax=max(responses*1.1)
    if ymin==None:
        ymin=min(indicies)*0.1
    plt.axis([min(indicies),max(indicies),ymin,ymax])
    plt.show()
    

def HarmonicResponseViewer(bw,harmonicResponse,title="",mode="COS_SIN"):
    #
    # okay I should make the matrix full of zeros and then set the values
    # as I go through the harmonicResponse

    #First we do regular bar charts plot

    #pprint(harmonicResponse)
    data=[[a,b,c,d] for (a,b),(c,d) in harmonicResponse]
    for mval in range(bw+1):
        HarmonicBarChart(mval,data,mchoiceMinusOne=False)
        HarmonicBarChart(mval,data,showPrime=True,mchoiceMinusOne=False)
    #Then the matrix plots

    
    mat_cos=zeros((bw,bw))
    mat_sin=mat_cos
    for (x,y),(c,s) in harmonicResponse:
        mat_cos[y][x]=c
        mat_sin[y][x]=s

    pylab.matshow(mat_cos)
    pylab.title(title+" cos")
    pylab.xlabel("n")
    pylab.ylabel("m")
    pylab.show()
    pylab.matshow(mat_sin)
    pylab.title(title+" sin")
    pylab.xlabel("n")
    pylab.ylabel("m")
    pylab.show()

    
def fieldplot(filename,
              BigMF_XYZSA=[True,True,True,True,True],
              SmallMF_XYZSA=[True,True,True,True,True],
              BigVP_XYZSA=[True,True,True,True,True],
              SmallVP_XYZSA=[True,True,True,True,True]):
    """We load a file like the one made by fieldfind"""

    def resultextractor(resultlist):
        """ returns (resultx,resulty,resultz,results) """
        resultx=[]
        resulty=[]
        resultz=[]
        results=[]
        for row in resultlist:
            tmpx=[p[0] for p in row]
            tmpy=[p[1] for p in row]
            tmpz=[p[2] for p in row]
            tmps=[p[3] for p in row]
            resultx.append(tmpx)
            resulty.append(tmpy)
            resultz.append(tmpz)
            results.append(tmps)
        return resultx,resulty,resultz,results


    f=file(filename,"r")
    q=cPickle.load(f)
    f.close()
    description=q[0]
    radialunits=q[1][0]
    axialunits=q[1][1]
    fieldunitsMF=q[1][2]
    vectpotunits=q[1][3]
    Big_rvectMF=q[2]
    Big_zvectMF=q[3]
    Big_rMF=q[4]

    Big_resultxMF,Big_resultyMF,Big_resultzMF,Big_resultsMF=resultextractor(Big_rMF)

    

 
    Small_rvectMF=q[5]
    Small_zvectMF=q[6]
    Small_rMF=q[7]
    Small_resultxMF,Small_resultyMF,Small_resultzMF,Small_resultsMF=resultextractor(Small_rMF)

    harmonicBW=q[8]
    RSphereList=q[9]
    BFieldonSpheres=q[10]
    BFieldSHTs=q[11]

    Big_rvectVP=q[12]
    Big_zvectVP=q[13]
    Big_rVP=q[14]

    Big_resultxVP,Big_resultyVP,Big_resultzVP,Big_resultsVP=resultextractor(Big_rVP)

    Small_rvectVP=q[15]
    Small_zvectVP=q[16]
    Small_rVP=q[17]

    Small_resultxVP,Small_resultyVP,Small_resultzVP,Small_resultsVP=resultextractor(Small_rVP)

    


    def fieldcontourplot(zvect,rvect,funcmap,title,numcontours=20,
                         linewidths=3,contourvalues=None,maxfield=None):
        #maxfield=0.00001
        if contourvalues == None and maxfield!=None:
            df=maxfield/numcontours
            contourvalues=arange(-maxfield,maxfield+df,df)
        if contourvalues==None and numcontours==None:
            raise Exception
        if contourvalues!=None:
            pylab.contour(zvect,rvect,funcmap,contourvalues,linewidths=linewidths)#fieldvals)
        else:
            pylab.contour(zvect,rvect,funcmap,numcontours,linewidths=linewidths)#fieldvals)
        pylab.axis('equal')
        pylab.title(title)
        pylab.ylabel("$r$ ("+radialunits+")")
        pylab.xlabel("$z$ ("+axialunits+")")
        pylab.show()

    def centralaxisplot(zvect,resultz,fieldunits):
        print "centralaxisplot"
        numR=len(resultz)

        if not is_odd(numR):
            print numR
            #raise Exception
        centralR=int((numR-1)/2)
        #print centralR
        numZ=len(resultz[0])
        centralZalongZ=resultz[centralR]
        if fieldunits=="teslas":
            centralZalongZ=1000.0*array(centralZalongZ)
            fieldunits="milliteslas"
        pylab.plot(zvect,centralZalongZ)
        pylab.title("$B_{z}$ field on $z$-axis")
        pylab.ylabel("$B_{z}$ ("+fieldunits+")")
        pylab.xlabel("$z$ ("+axialunits+")")
        pylab.show()

    def radialaxisplot(rvect,  # a list of the radial positions
                       resultx,
                       resulty,
                       resultz,
                       fieldunits):
        print "radialaxisplot"
        #print "rvect:"
        #pprint(rvect)
        
        numZ=len(resultz[0])
        numR=len(resultz)
        print "resultsz:(",numR,",",numZ,")"
        if not is_odd(numZ):
            print numZ
        centralZ=int((numZ-1)/2)
        centralR=int((numR-1)/2)
        centralRalongR=[resultz[i][centralZ] for i in range(numZ)]
        if fieldunits=="teslas":
            centralRalongR=1000.0*array(centralRalongR)
            fieldunits="milliteslas"
        pylab.plot(rvect,centralRalongR)
        pylab.title("$B_{z}$ field on $r$-axis")
        pylab.ylabel("$B_{z}$ ("+fieldunits+")")
        pylab.xlabel("$r$ ("+radialunits+")")
        pylab.show()
        
        centralRalongR=[resultx[i][centralZ] for i in range(numZ)]
        if fieldunits=="teslas":
            centralRalongR=1000.0*array(centralRalongR)
            fieldunits="milliteslas"
        pylab.plot(rvect,centralRalongR)
        pylab.title("$B_{x}$ field on $x$-axis")
        pylab.ylabel("$B_{x}$ ("+fieldunits+")")
        pylab.xlabel("$x$ ("+radialunits+")")
        pylab.show()

        centralRalongR=[resulty[i][centralZ] for i in range(numZ)]
        if fieldunits=="teslas":
            centralRalongR=1000.0*array(centralRalongR)
            fieldunits="milliteslas"
        pylab.plot(rvect,centralRalongR)
        pylab.title("$B_{y}$ field on $x$-axis")
        pylab.ylabel("$B_{y}$ ("+fieldunits+")")
        pylab.xlabel("$x$ ("+radialunits+")")
        pylab.show()

        
        
        
        


    print "BIGRESULTY [10,10]:",Big_resultyMF[0][0]
    #Big_XYZSA=Small_XYZSA=[False for i in range(5)]

    # Axis plots
    if BigMF_XYZSA[4]:
        radialaxisplot(Big_zvectMF,
                       Big_resultxMF,Big_resultyMF,Big_resultzMF,
                       fieldunitsMF)
        centralaxisplot(Big_rvectMF,Big_resultzMF,fieldunitsMF)
    if SmallMF_XYZSA[4]:
        radialaxisplot(Small_rvectMF,
                       Small_resultxMF,Small_resultyMF,Small_resultzMF,
                       fieldunitsMF)
        centralaxisplot(Small_zvectMF,Small_resultzMF,fieldunitsMF)


    # Magnetic Field
    if BigMF_XYZSA[0]:
        fieldcontourplot(Big_zvectMF,Big_rvectMF,Big_resultxMF,"$B_{x}$ field")
    if SmallMF_XYZSA[0]:
        fieldcontourplot(Small_zvectMF,Small_rvectMF,Small_resultxMF,"$B_{x}$ field")
    if BigMF_XYZSA[1]:
        fieldcontourplot(Big_zvectMF,Big_rvectMF,Big_resultyMF,"$B_{y}$ field")
    if SmallMF_XYZSA[1]:
        fieldcontourplot(Small_zvectMF,Small_rvectMF,Small_resultyMF,"$B_{y}$ field")
        
    if BigMF_XYZSA[2]:
        fieldcontourplot(Big_zvectMF,Big_rvectMF,Big_resultzMF,"$B_{z}$ field")
    if SmallMF_XYZSA[2]:
        fieldcontourplot(Small_zvectMF,Small_rvectMF,Small_resultzMF,"$B_{z}$ field")
    if BigMF_XYZSA[3]:
        fieldcontourplot(Big_zvectMF,Big_rvectMF,Big_resultsMF,"$B^{2}$")
    if SmallMF_XYZSA[3]:
        fieldcontourplot(Small_zvectMF,Small_rvectMF,Small_resultsMF,"$B^{2}$")
    
    


    # Vector Potential
    if BigVP_XYZSA[0]:
        fieldcontourplot(Big_zvectVP,Big_rvectVP,Big_resultxVP,"$A_{x}$ field")
    if SmallVP_XYZSA[0]:
        fieldcontourplot(Small_zvectVP,Small_rvectVP,Small_resultxVP,"$A_{x}$ field")
    if BigVP_XYZSA[1]:
        fieldcontourplot(Big_zvectVP,Big_rvectVP,Big_resultyVP,"$A_{y}$ field")
    if SmallVP_XYZSA[1]:
        fieldcontourplot(Small_zvectVP,Small_rvectVP,Small_resultyVP,"$A_{y}$ field")
        
    if BigVP_XYZSA[2]:
        fieldcontourplot(Big_zvectVP,Big_rvectVP,Big_resultzVP,"$A_{z}$ field")
    if SmallVP_XYZSA[2]:
        fieldcontourplot(Small_zvectVP,Small_rvectVP,Small_resultzVP,"$A_{z}$ field")
    if BigVP_XYZSA[3]:
        fieldcontourplot(Big_zvectVP,Big_rvectVP,Big_resultsVP,"$A^{2}$")
    if SmallVP_XYZSA[3]:
        fieldcontourplot(Small_zvectVP,Small_rvectVP,Small_resultsVP,"$A^{2}$")
    
    #if BigVP_XYZSA[4]:
    #    centralaxisplot(Big_zvectVP,Big_resultzVP,vectpotunits)
    #if SmallVP_XYZSA[4]:
    #    centralaxisplot(Small_zvectVP,Small_resultzVP,vectpotunits)


    # Spherical Harmonic Plot
    for RSphere,HRS in zip(RSphereList,BFieldSHTs):
        # Note that we assume the Harmonic Response has been computed with Rref=1.0
        HarmonicResponseViewer(harmonicBW,HRS[2],title="$B_{z}$ "+RSphere.__str__())
        HarmonicResponseViewer(harmonicBW,HRS[0],title="$B_{x}$ "+RSphere.__str__())
        HarmonicResponseViewer(harmonicBW,HRS[1],title="$B_{y}$ "+RSphere.__str__())
        HarmonicResponseViewer(harmonicBW,HRS[3],title="$B^{2}$ "+RSphere.__str__())
