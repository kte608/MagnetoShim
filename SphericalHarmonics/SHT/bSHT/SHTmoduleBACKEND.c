#include "SHTmodule.h"


/*static PyObject*
SphericalHarmonicTransform(PyObject *self,PyObject *args)*/


void denormalize(double *rresult,double* iresult,int bw);

void SHT_BACKEND(long int bw,double *rdata,double *idata,
		 double **rresult,double **iresult)
{

  //long int bw;
  //int i, size;
  int size;
  int cutoff ;
  int rank, howmany_rank ;
  double *rcoeffs, *icoeffs;
  double *workspace, *weights;
  fftw_plan dctPlan;
  fftw_plan fftPlan;
  fftw_iodim dims[1], howmany_dims[1];


  //PyObject* result;
  //PyObject* tempr;
  //PyObject* tempi;
  // PyObject* tempobject;


  /*//Take in arguments
  if(!PyArg_ParseTuple(args,"O",
		       c))
  {return NULL;}
  Py_INCREF(c);
  b=PyTuple_GetItem(c,0);
  Py_INCREF(b);
  bw=PyInt_AsLong(PyTuple_GetItem(b,0));
  realdata=PyTuple_GetItem(b,1);
  Py_INCREF(realdata);
  imagdata=PySequence_GetItem(b,2);
  Py_INCREF(imagdata);*/


  /*//convert the python lists to double arrays
  rdata=listtoarray(realdata);
  idata=listtoarray(imagdata);*/

    /*** ASSUMING WILL SEMINAIVE ALL ORDERS ***/
  cutoff = bw ;
  size = 2*bw;

  if((rdata==NULL)||(idata==NULL)){printf("Bad Data in SHTBACKEND\n");}
  if(rdata==NULL){printf("Bad rdata\n");}
  if(idata==NULL){printf("Bad idata\n");}
  //printf("About to allocate Memory\n");
  rcoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  icoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  *rresult = (double *) malloc(sizeof(double) * (bw * bw));
  *iresult = (double *) malloc(sizeof(double) * (bw * bw));
  workspace = (double *) malloc(sizeof(double) * 
				((10 * (bw*bw)) + 
				 (24 * bw)));

   // make array for weights 
  weights = (double *) malloc(sizeof(double) * 4 * bw);
  //printf("Done Allocating Memory\n");


  
  //  At this point, check to see if all the memory has been
  //  allocated. If it has not, there's no point in going further.
  //printf("About to allocate Memory\n");
  if ( (rdata == NULL) || (idata == NULL) ||
       (*rresult == NULL) || (*iresult == NULL) ||
       (rcoeffs == NULL) || (icoeffs == NULL) ||
       (workspace == NULL)||(weights == NULL))
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }



  /*
    construct fftw plans
  */

  /* make DCT plans -> note that I will be using the GURU
     interface to execute these plans within the routines*/

  /* forward DCT */
  dctPlan = fftw_plan_r2r_1d( 2*bw, weights, rdata,
			      FFTW_REDFT10, FFTW_ESTIMATE ) ;  
  /*
    fftw "preamble" ;
    note that this plan places the output in a transposed array
  */
  rank = 1 ;
  dims[0].n = 2*bw ;
  dims[0].is = 1 ;
  dims[0].os = 2*bw ;
  howmany_rank = 1 ;
  howmany_dims[0].n = 2*bw ;
  howmany_dims[0].is = 2*bw ;
  howmany_dims[0].os = 1 ;
  /* forward fft */
  fftPlan = fftw_plan_guru_split_dft( rank, dims,
				      howmany_rank, howmany_dims,
				      rdata, idata,
				      workspace, workspace+(4*bw*bw),
				      FFTW_ESTIMATE ); 

  /* now make the weights */
  makeweights( bw, weights );


   /* now do the forward spherical transform */
   // tstart = csecond();

    FST_semi_fly(rdata, idata,
		 *rresult, *iresult,
		 bw,
		 workspace,
		 0,
		 cutoff,
		 &dctPlan,
		 &fftPlan,
		 weights ) ;
  /* destroy fftw plans */
  fftw_destroy_plan( fftPlan );
  fftw_destroy_plan( dctPlan );

  // In this section the results are modified so that they represent
  // the un-normalized spherical harmonics.


  //printf("I am here\n");
  denormalize(*rresult,*iresult,bw);
  //printf("Done Denormalization\n");
 

 /* free memory */
  free( weights );
  free(workspace);
  //free(iresult);
  // free(rresult);
  //free(idata);
  //free(rdata);
  free(icoeffs);
  free(rcoeffs);

  /*// construct the result list
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
  return Py_BuildValue("O",result);*/
  return;
}



void denormalize(double *rresult,double* iresult,int bw)
{
 // In this section the results are modified so that they represent
  // the un-normalized spherical harmonics.
  int i,j,index;
    double temp,coeff;
    //printf("Started Denormalization\n");
  for(i=0;i<bw;i++)
  { 
      coeff=harmcoeff(i,0);
      index=seanindex(0,i,bw);
      temp=rresult[index];
      temp=temp*coeff;

      rresult[index]=temp;
      temp=iresult[index];
      temp=temp*coeff;
      iresult[index]=temp;
      for(j=1;j<=i;j++)
      {
	  coeff=harmcoeff(i,j);
	  index=seanindex(j,i,bw);
	  temp=rresult[index];
	  temp=temp*coeff;
	  rresult[index]=temp;

	  temp=iresult[index];
	  temp=temp*coeff;
	  iresult[index]=temp; 

	  index=seanindex(-j,i,bw);
	  temp=rresult[index];
	  temp=temp*coeff;
	  rresult[index]=temp;

	  temp=iresult[index];
	  temp=temp*coeff;
	  iresult[index]=temp;
      }
  }


}

/*
static PyObject*
InvSphericalHarmonicTransform(PyObject *self,PyObject *args)
{
  FILE *errorsfp;
  int i, j, bw, size, loops;
  int l, m, dummy, cutoff ;
  int rank, howmany_rank ;
  double *rcoeffs, *icoeffs, *rdata, *idata, *rresult, *iresult;
  double *workspace, *weights;
  double dumx, dumy ;
  double *relerror, *curmax, granderror, grandrelerror;
  double realtmp, imagtmp,origmag, tmpmag;
  double ave_error, ave_relerror, stddev_error, stddev_relerror;
  double total_time, for_time, inv_time;
  double tstart, tstop;
  time_t seed;
  fftw_plan dctPlan, idctPlan ;
  fftw_plan fftPlan, ifftPlan ;
  fftw_iodim dims[1], howmany_dims[1];


  // ASSUMING WILL SEMINAIVE ALL ORDERS **
  cutoff = bw ;

  size = 2*bw;
  total_time = 0.0;
  for_time = 0.0;
  inv_time = 0.0;
  granderror = 0.0;
  grandrelerror = 0.0; 


 
  //  allocate memory
 

  rcoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  icoeffs = (double *) malloc(sizeof(double) * (bw * bw));
  rdata = (double *) malloc(sizeof(double) * (size * size));
  idata = (double *) malloc(sizeof(double) * (size * size));
  rresult = (double *) malloc(sizeof(double) * (bw * bw));
  iresult = (double *) malloc(sizeof(double) * (bw * bw));
  workspace = (double *) malloc(sizeof(double) * 
				((10 * (bw*bw)) + 
				 (24 * bw)));

  //space for errors 
  relerror = (double *) malloc(sizeof(double) * loops);
  curmax = (double *) malloc(sizeof(double) * loops);

  // make array for weights /
  // weights = (double *) malloc(sizeof(double) * 4 * bw);

  //
  //  At this point, check to see if all the memory has been
  //  allocated. If it has not, there's no point in going further.
 
  if ( (rdata == NULL) || (idata == NULL) ||
       (rresult == NULL) || (iresult == NULL) ||
       (rcoeffs == NULL) || (icoeffs == NULL) ||
       (workspace == NULL) || (weights == NULL) )
    {
      perror("Error in allocating memory");
      exit( 1 ) ;
    }

  
  //  construct fftw plans
  

  // make DCT plans -> note that I will be using the GURU
 //  interface to execute these plans within the routines

  // forward DCT 
  dctPlan = fftw_plan_r2r_1d( 2*bw, weights, rdata,
			      FFTW_REDFT10, FFTW_ESTIMATE ) ;
      
  // inverse DCT 
  idctPlan = fftw_plan_r2r_1d( 2*bw, weights, rdata,
			       FFTW_REDFT01, FFTW_ESTIMATE );

  
  //  fftw "preamble" ;
  //  note that this plan places the output in a transposed array
  
  rank = 1 ;
  dims[0].n = 2*bw ;
  dims[0].is = 1 ;
  dims[0].os = 2*bw ;
  howmany_rank = 1 ;
  howmany_dims[0].n = 2*bw ;
  howmany_dims[0].is = 2*bw ;
  howmany_dims[0].os = 1 ;
  
  // forward fft 
  fftPlan = fftw_plan_guru_split_dft( rank, dims,
				      howmany_rank, howmany_dims,
				      rdata, idata,
				      workspace, workspace+(4*bw*bw),
				      FFTW_ESTIMATE );
  
  /
  //now plan for inverse fft - note that this plans assumes
  //  that I'm working with a transposed array, e.g. the inputs
  // for a length 2*bw transform are placed every 2*bw apart,
  // the output will be consecutive entries in the array
  
  rank = 1 ;
  dims[0].n = 2*bw ;
  dims[0].is = 2*bw ;
  dims[0].os = 1 ;
  howmany_rank = 1 ;
  howmany_dims[0].n = 2*bw ;
  howmany_dims[0].is = 1 ;
  howmany_dims[0].os = 2*bw ;

  // inverse fft 
  ifftPlan = fftw_plan_guru_split_dft( rank, dims,
				       howmany_rank, howmany_dims,
				       rdata, idata,
				       workspace, workspace+(4*bw*bw),
				       FFTW_ESTIMATE );

  // now make the weights 
  makeweights( bw, weights );


   // now do the forward spherical transform 
    tstart = csecond();

    FST_semi_fly(rdata, idata,
		 rresult, iresult,
		 bw,
		 workspace,
		 1,
		 cutoff,
		 &dctPlan,
		 &fftPlan,
		 weights ) ;

    tstop = csecond();    
    for_time += (tstop - tstart);


    // do the inverse spherical transform 
    tstart = csecond();    

    InvFST_semi_fly(rcoeffs,icoeffs,
		    rdata, idata,
		    bw,
		    workspace,
		    1,
		    cutoff,
		    &idctPlan,
		    &ifftPlan );

    tstop = csecond();
    inv_time += (tstop - tstart);


  // destroy fftw plans 
  fftw_destroy_plan( ifftPlan );
  fftw_destroy_plan( fftPlan );
  fftw_destroy_plan( idctPlan );
  fftw_destroy_plan( dctPlan );

  // free memory 
  free( weights );
  free(curmax);
  free(relerror);
  free(workspace);
  free(iresult);
  free(rresult);
  free(idata);
  free(rdata);
  free(icoeffs);
  free(rcoeffs);

  return 0 ;



}
*/
