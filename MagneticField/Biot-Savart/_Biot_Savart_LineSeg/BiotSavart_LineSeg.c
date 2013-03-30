#include <math.h>
#include <Python.h>
#include <pthread.h>    //for multithreading
#include <semaphore.h>  //for multithreading
#include "Vector.h"
#include "LineSeg.h"
const double mu_0 = 4*M_PI*1e-7; // N/A^2
static int numthreads=1;

Vector Biot_Savart_SingleLineSeg(LineSeg seg,Vector fieldpoint)
{

    Vector result,r,Seg_Cross_r;
    Vector SegMidPoint,SegVector,SegDirection;
    Vector BDirection;
    double SegLength,x,z,constval;
    double Bmag;


    SegMidPoint=linesegCenter(seg);
    SegLength=linesegLength(seg);
    SegVector=linesegVector(seg);
    SegDirection=vectorUnitVector(SegVector);


    r=linesegVector(linesegGenerate(SegMidPoint,fieldpoint));
    Seg_Cross_r=vectorCross(SegVector,r);
    BDirection=vectorUnitVector(Seg_Cross_r);
    z=vectorMagnitude(Seg_Cross_r)/SegLength;
    constval=mu_0/(4*M_PI*z);
    x=vectorDot(SegDirection,r);
    if(z == 0)
    {
	Bmag=0;
    }
    else
    {
	Bmag=constval*((SegLength-2*x)/sqrt((SegLength-2*x)*(SegLength-2*x)+4*z*z)+(SegLength+2*x)/sqrt((SegLength+2*x)*(SegLength+2*x)+4*z*z));
    }
    //printf("Bmag %g\n",Bmag);
    result=vectorScale(BDirection,Bmag);
    //printf("result.x,result.y,result.z %g,%g,%g\n",result.x,result.y,result.z);
    return result;
}

typedef struct
{
  int start;
  int end;
  LineSeg *segs;
  int numsegs;
  Vector* fieldpoints;
  Vector* Bfield;
} WorkerArg; // Nice for threads

void *MagFieldLineSegArray_BACKENDWORKER(void *arg)
{
  WorkerArg Arg;
  int fpFirst,fpLast;
  int fp,seg;
  
  Vector zerovect;
  Vector tmpvectA;
  Vector tmpvectB;

  double percent;
  int percent_denom;
  zerovect.coords[0]=0.0;
  zerovect.coords[1]=0.0;
  zerovect.coords[2]=0.0;
  
  Arg=*((WorkerArg*)(arg));
  fpFirst=Arg.start;
  fpLast=Arg.end;

  for(fp=fpFirst;fp<=fpLast;fp++)
    {
      percent_denom=(fpLast-fpFirst);
      if(percent_denom==0)
	{percent=100.0;}
      else
	{percent=100.0*(fp-fpFirst)/percent_denom;}
      printf("Fieldpoint %d in %d to %d. %g %% complete.\n",fp,fpFirst,fpLast,percent);

      
      tmpvectA=zerovect;
      for(seg=0;seg<Arg.numsegs;seg++)
	{
	  
	    tmpvectB=vectorAdd(tmpvectA,Biot_Savart_SingleLineSeg(Arg.segs[seg],Arg.fieldpoints[fp]));  
	    tmpvectA=tmpvectB;
	
	}
      Arg.Bfield[fp]=tmpvectB;
    }
  return NULL;
}


Vector *MagFieldLineSegArray_BACKEND(LineSeg *segs,int numsegs,
				     Vector* fieldpoints,int numpoints)
{
    int i;
    Vector* Bfield;
    WorkerArg *args;
    pthread_t *tid;
    int err;
    int fieldPointsPerThread;
    Bfield=malloc(sizeof(Vector)*numpoints);
    if(Bfield==NULL){printf("Bfield allocation failure\n"); exit(-1);}

    tid=malloc(sizeof(pthread_t)*numthreads);
    if(tid==NULL){printf("tid allocation failed\n"); exit(-1);}
    args=malloc(sizeof(WorkerArg)*numthreads);
    if(args==NULL){printf("arg allocation failed\n"); exit(-1);}

    // New stuff
    fieldPointsPerThread=numpoints/numthreads;
    printf("Starting threads. fieldpoints:%d fieldPointsPerThread:%d\n",numpoints,fieldPointsPerThread);
    for(i=0;i<numthreads;i++)
      {
	args[i].start=i*fieldPointsPerThread;
      if(i!=(numthreads-1))
	{
	  args[i].end=(i+1)*fieldPointsPerThread-1;
	}
      else
	{
	  // We need to take care of the remainder too..
	  args[i].end=numpoints-1;
	}
      args[i].segs=segs;
      args[i].Bfield=Bfield;
      args[i].fieldpoints=fieldpoints;
      args[i].numsegs=numsegs;
      err=pthread_create(&(tid[i]),NULL,
			 &MagFieldLineSegArray_BACKENDWORKER,
			 (void*)&(args[i]));
      if(err != 0)
	{printf("Failed to create thread\n"); exit(-1);}
    }

  // Joins.
  for(i=0;i<numthreads;i++)
    {
      pthread_join(tid[i],NULL);
    }
  
  free(tid);
  free(args);



  return Bfield;
}





double* listtoarray(PyObject* o,int *size)
{
    int length;
    double* result;
    double tempval;
    PyObject* tempObject;
    int i;
    //printf("Started List to Array\n");
    length = PySequence_Size(o);
    //printf("%d\n",length);
    result=malloc(length*sizeof(double));
    for(i=0;i<length;i++)
    {
      //printf("%d\n",i);
	tempObject=PySequence_GetItem(o,i);
	tempval=PyFloat_AsDouble(tempObject);
	//printf("%g\n",tempval);
	result[i]=tempval;
    }
    //printf("The length is: %d\n",length);
    *size=length;
    return result;
}



static PyObject* MagFieldLineSegArray(PyObject *self,PyObject *args)
{
// As an argument I am expecting first a list of LineSegs of the form
// [[[xs,ys,zs],[xe,ye,ze]],...] which describe the coil
// I am then expecting a list of Vectors of the form [[x,y,z],...]
// which describe the fieldpoints
// The returned value will be a list of field vectors of the form [[x,y,z,mag2],...]
    //double x,y,z;
    Vector A;
    Vector B;
    int i;//,tmpint;
    PyObject *SegList;
    PyObject *tmpSeg;
    PyObject *start,*end;
//    PyObject *p_Seg;
    LineSeg tempseg;
    int numsegs;
    LineSeg* segs;
    
    PyObject *FieldPoints;
    PyObject *p_point;
    int numpoints;
    Vector *points;

    Vector *Bfield;
    PyObject *result; // the magnetic field in [[x,y,z,mag2],...] format
    PyObject *elem;   // [x,y,z,mag2]
    // Find the SegList and FieldPoints lists in Python Format

    if (!PyArg_ParseTuple(args, "OO", &SegList,&FieldPoints))
        return NULL;

    // Convert SegList from python format to segs in C format
    numsegs=PySequence_Size(SegList);
    segs=malloc(numsegs*sizeof(LineSeg));
    //printf("numsegs %d\n",numsegs);
    for(i=0;i<numsegs;i++)
    {
	tmpSeg=PySequence_GetItem(SegList,i);
	if(PySequence_Size(tmpSeg) != 2)
	{printf("Invalid Seg Format \n"); return NULL;}
	start=PySequence_GetItem(tmpSeg,0);
	A.coords[0]=PyFloat_AsDouble(PySequence_GetItem(start,0));
	A.coords[1]=PyFloat_AsDouble(PySequence_GetItem(start,1));
	A.coords[2]=PyFloat_AsDouble(PySequence_GetItem(start,2));
	end=PySequence_GetItem(tmpSeg,1);
	B.coords[0]=PyFloat_AsDouble(PySequence_GetItem(end,0));
	B.coords[1]=PyFloat_AsDouble(PySequence_GetItem(end,1));
	B.coords[2]=PyFloat_AsDouble(PySequence_GetItem(end,2));
	tempseg = linesegGenerate(A,B);
	segs[i]=tempseg;
    }

    // Convert FieldPoints from python format to points in C format
    numpoints=PySequence_Size(FieldPoints);
    points=malloc(numpoints*sizeof(Vector));
    //printf("numpoints %d\n",numpoints);
    for(i=0;i<numpoints;i++)
    {
	p_point=PySequence_GetItem(FieldPoints,i);
	if(PySequence_Size(p_point) != 3)
	{printf("Invalid fieldpoint format\n"); return NULL;}
	A.coords[0]=PyFloat_AsDouble(PySequence_GetItem(p_point,0));
	A.coords[1]=PyFloat_AsDouble(PySequence_GetItem(p_point,1));
	A.coords[2]=PyFloat_AsDouble(PySequence_GetItem(p_point,2));
	points[i]=A;
    }

    //printf("I parsed the arguments\n");
    // do the work
    //printf("About to Start the Magnetic Field Kernel\n");
    Bfield=MagFieldLineSegArray_BACKEND(segs,numsegs,
					points,numpoints);
    //printf("Finished computing the Field, returning result to Python\n");

    // Convert the Bfield from C format to Python format and add the mag2 
    result=PyList_New(numpoints);

    for(i=0;i<numpoints;i++)
    {
	B=Bfield[i];
	//printf("Bfield[i].x %g\n",B.x);
	elem=PyList_New(4);
	PyList_SetItem(elem,0,PyFloat_FromDouble(B.coords[0]));
	PyList_SetItem(elem,1,PyFloat_FromDouble(B.coords[1]));
	PyList_SetItem(elem,2,PyFloat_FromDouble(B.coords[2]));
	PyList_SetItem(elem,3,PyFloat_FromDouble(B.coords[0]*B.coords[0]+B.coords[1]*B.coords[1]+B.coords[2]*B.coords[2]));
	PyList_SetItem(result,i,elem);
    }
    // Return
    return result;
}



static PyObject* SetNumThreads(PyObject *self,PyObject *args)
{
  PyObject *result;
  
  if (!PyArg_ParseTuple(args, "i", &numthreads))
    return NULL;
  
  result=Py_BuildValue("");
  return result;
}

static PyObject *BiotSavartLineSeg_Error;



static PyMethodDef BiotSavartLineSegMethods[]={
  {"MagFieldLineSegArray",MagFieldLineSegArray,METH_VARARGS,
   "LineSegs,FieldPoints"},
  {"SetNumThreads",SetNumThreads,METH_VARARGS,"numthreads"},
  {NULL,NULL,0,NULL}/*Sentinel*/
};

PyMODINIT_FUNC
initBiotSavartLineSeg(void)
{
  PyObject *m;

  m=Py_InitModule("BiotSavartLineSeg",BiotSavartLineSegMethods);

  BiotSavartLineSeg_Error = PyErr_NewException("BiotSavartLineSeg.error",NULL,NULL);

  Py_INCREF(BiotSavartLineSeg_Error);
  PyModule_AddObject(m,"error",BiotSavartLineSeg_Error);
}
