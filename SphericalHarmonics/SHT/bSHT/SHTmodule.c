#include "SHTmodule.h"

static PyObject*
SHT_Base_Level_Test(PyObject *self,PyObject *args)
{
  // This function runs a couple of simple tests on the backend and prints
  // the results to the screen.
  int i;

  double *rresult,*iresult; 
  long int bw=2;
  double rdata[]={1.0,1.0,1.0,1.0,
		  1.0,1.0,1.0,1.0,
		  1.0,1.0,1.0,1.0,
		  1.0,1.0,1.0,1.0};
  double idata[]={0.0,0.0,0.0,0.0,
		  0.0,0.0,0.0,0.0,
		  0.0,0.0,0.0,0.0,
		  0.0,0.0,0.0,0.0};
  SHT_BACKEND(bw,rdata,idata,&rresult,&iresult);

  for(i=0;i<4;i++)
    {
      printf("Real: %g\t\t\t\tImag: %g\n",rresult[i],iresult[i]);
      }


  return Py_BuildValue("");
}


static PyObject*
SphericalHarmonicTransform(PyObject *self,PyObject *args)

{
  PyObject* b;
  PyObject* c;

  PyObject* realdata;
  PyObject* imagdata;
  double *rdata,*idata;
  double *rresult,*iresult;
  long int bw;
  int i;
  int datasize;
  PyObject* result;
  PyObject* tempr;
  PyObject* tempi;
  PyObject* tempobject;
//Take in arguments
  Py_INCREF(args);
  //Py_INCREF(self);
  //printf("I am in the Spherical Harmonic transform Code\n");
  //printf("%s\n",PyString_AsString(PyObject_Str(args)));

  //if(!PyArg_ParseTuple(args,"O",
//		       c))
  // {return NULL;}
  //c=args;
  //Py_INCREF(c);
  //printf("%s\n",PyString_AsString(PyObject_Str(c)));
  //printf("Now I am here\n");

  b=PyTuple_GetItem(args,0); 
  //Py_INCREF(b);
  //printf("%s\n",PyString_AsString(PyObject_Str(b)));
  bw=PyInt_AsLong(PySequence_GetItem(b,0));
  realdata=PySequence_GetItem(b,1);
  //Py_INCREF(realdata);
  imagdata=PySequence_GetItem(b,2);
  //Py_INCREF(imagdata);
  //printf("%s\n",PyString_AsString(PyObject_Str(realdata)));
  //printf("%s\n",PyString_AsString(PyObject_Str(imagdata)));
  //printf("I am here\n");
  //convert the python lists to double arrays
  rdata=listtoarray(realdata,&datasize);
  //printf("Datasize: %d\n",datasize);
  //for(i=0;i<datasize;i++){printf("%g\n",rdata[i]);}
  idata=listtoarray(imagdata,&datasize);
  //for(i=0;i<datasize;i++){printf("%g\n",idata[i]);}
  if((rdata==NULL)||(idata==NULL)){printf("Bad Data\n");}
  //printf("About to do the backend stuff\n");
  SHT_BACKEND(bw,rdata,idata,&rresult,&iresult);
  //printf("About to construct result list\n");

 // construct the result list
  result=PyList_New(bw*bw);
  for(i=0;i<bw*bw;i++)
  {
      tempr=PyFloat_FromDouble(rresult[i]);
      tempi=PyFloat_FromDouble(iresult[i]);
      tempobject=PyTuple_New(2);
      PyTuple_SetItem(tempobject,0,tempr);
      PyTuple_SetItem(tempobject,1,tempi);
      PyList_SetItem(result,i,tempobject);
  }
  //Py_DECREF(b);
  //Py_DECREF(c);
  return Py_BuildValue("O",result);

  
}



double factorial(int n)
{
    double result=1;
    int i;
    for(i=n;i>1;i--)
    {
	result=result*i;
    }
    return result;
}

double harmcoeff(int n,int m)
{
    double temp1,temp2,temp3;
    temp1=factorial(n-m);
    temp2=factorial(n+m);
    temp3=temp1/temp2;
    //printf("%g %g %g\t\t",temp1,temp2,temp3);
    return sqrt(((2*n+1)*temp3)/((4*M_PI)));
}

/**
 * \brief
 * converts a python list into an array
 * \param o
 * the python object which should be a list
 * \param *size
 * this is loaded with the size of the resulting array
 */
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





static PyObject *bSHT_Error;



static PyMethodDef bSHTMethods[]={
  {"SHT",SphericalHarmonicTransform,METH_VARARGS,
   "(bw,realdata,imagdata)"},
  {"SHT_Base_Level_Test",SHT_Base_Level_Test,METH_VARARGS,
   "This does not take any arguments but it does use printf to print to the screen"},
 
 
  {NULL,NULL,0,NULL}/*Sentinel*/
};

PyMODINIT_FUNC
initbSHT(void)
{
  PyObject *m;
  //printf("In SHTmodule.c .... Initializing bSHT\n");
  m=Py_InitModule("bSHT",bSHTMethods);

  bSHT_Error = PyErr_NewException("bSHT.error",NULL,NULL);

  Py_INCREF(bSHT_Error);
  PyModule_AddObject(m,"error",bSHT_Error);
}
