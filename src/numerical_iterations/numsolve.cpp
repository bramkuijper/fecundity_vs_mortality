#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "auxiliary.h"

using namespace std;

string filename("iter_fecmort");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

struct rparams
{
    // VARS
};


double bound0(double val)
{
    if (val < 0)
    {
        return(0.0001);
    }

    return(val);
}

double bound01(double val)
{
    val = val < 0 ? 0.0001 : val > 1.0 ? 0.9999 : val;

    return(val);
}

// recursions of all the patch frequencies
int psys_recur(void *params, gsl_vector *f)
{
    //VARFUNC
    
    
    // PATCHRECUR
    
    return GSL_SUCCESS;
}


// reproductive values
int rvsys_recur(void *params, gsl_vector *f)
{
    //VARFUNC
    
    
    // REPVALRECUR
    
    return GSL_SUCCESS;
}

// relatedness coefficients
int relcoeff_recur(void *params, gsl_vector *f)
{
    //VARFUNC
   
    // RELRECUR 
    
    return GSL_SUCCESS;
}

void selgrads(void *params, gsl_vector *f)
{
    // VARFUNC

    // SELGRADS

}

void write_params(void *params)
{
    // VARFUNC


    // WRITEPARAMS
}


void write_data(void *params, int time)
{
    // VARFUNC


    if (time < 0)
    {
        // HEADERWRT
    }

    // WRITEDATA
}

int main (int argc, char **argv)
{
    long int max_iter = atoi(argv[1]);

    // initialize the vectors that contain the variables
    // functions solve for
    gsl_vector *x_p = gsl_vector_alloc(6);
    gsl_vector *x_rv = gsl_vector_alloc(8);
    gsl_vector *x_rel = gsl_vector_alloc(6);
    gsl_vector *x_selgrad = gsl_vector_alloc(3);

    // initialize the struct with parameters
    struct rparams paramstruct; 
    
    // ARGINIT


    write_data(&paramstruct,-1);
    // iterate
    for (int iter = 0; iter < max_iter ; ++iter)
    {
        psys_recur(&paramstruct, x_p);

        
        paramstruct.f10 = gsl_vector_get(x_p,0);
        paramstruct.f11 = gsl_vector_get(x_p,1);
        paramstruct.f12 = gsl_vector_get(x_p,2);
        paramstruct.f20 = gsl_vector_get(x_p,3);
        paramstruct.f21 = gsl_vector_get(x_p,4);
        paramstruct.f22 = gsl_vector_get(x_p,5);
        // reproductive values
        rvsys_recur(&paramstruct, x_rv);


        paramstruct.v1a0 = gsl_vector_get(x_rv,0);
        paramstruct.v1a1 = gsl_vector_get(x_rv,1);
        paramstruct.v1m0 = gsl_vector_get(x_rv,2);
        paramstruct.v1m1 = gsl_vector_get(x_rv,3);
        paramstruct.v2a0 = gsl_vector_get(x_rv,4);
        paramstruct.v2a1 = gsl_vector_get(x_rv,5);
        paramstruct.v2m0 = gsl_vector_get(x_rv,6);
        paramstruct.v2m1 = gsl_vector_get(x_rv,7);

        // relatedness
        relcoeff_recur(&paramstruct, x_rel);

        paramstruct.raa12 = gsl_vector_get(x_rel, 0);
        paramstruct.raa22 = gsl_vector_get(x_rel, 1);
        paramstruct.ram11 = gsl_vector_get(x_rel, 2);
        paramstruct.ram21 = gsl_vector_get(x_rel, 3);
        paramstruct.rmm10 = gsl_vector_get(x_rel, 4);
        paramstruct.rmm20 = gsl_vector_get(x_rel, 5);
        
        // selection gradients
        selgrads(&paramstruct, x_selgrad);

        bool condition_p1 = (fabs(paramstruct.p1 - gsl_vector_get(x_selgrad, 0)) < 1e-10 || paramstruct.p1 >= 0.999) || paramstruct.p1 <= 0.001;
        bool condition_p2 = (fabs(paramstruct.p2 - gsl_vector_get(x_selgrad, 1)) < 1e-10 || paramstruct.p2 >= 0.999) || paramstruct.p2 <= 0.001;

//        if (condition_p1 && condition_p2 &&
//            fabs(paramstruct.d - gsl_vector_get(x_selgrad, 2)) < 1e-07)
        if (condition_p1 && condition_p2)
        {
            paramstruct.p1 = gsl_vector_get(x_selgrad, 0);
            paramstruct.p2 = gsl_vector_get(x_selgrad, 1);
            paramstruct.d = paramstruct.d;//gsl_vector_get(x_selgrad, 2);

            write_data(&paramstruct,iter);
            break;
        }
        paramstruct.p1 = gsl_vector_get(x_selgrad, 0);
        paramstruct.p2 = gsl_vector_get(x_selgrad, 1);
        paramstruct.d = paramstruct.d;//gsl_vector_get(x_selgrad, 2);

        if ((isnan(paramstruct.p1) != 0 || isnan(paramstruct.p2) != 0) || isnan(paramstruct.d) != 0)
        {
            write_data(&paramstruct,iter);
            break;
        }

        if (iter % 10 == 0)
        {
            write_data(&paramstruct,iter);
        }
    }

    write_params(&paramstruct);


    gsl_vector_free(x_rel);
    gsl_vector_free(x_p);
    gsl_vector_free(x_rv);
    gsl_vector_free(x_selgrad);
}
