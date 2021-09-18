#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "c3/array.h"
#include "c3/lib_clinalg.h"
#include "c3/lib_funcs.h"

#include "c3/c3_interface.h"

#include <stdio.h>
#include <complex.h>
#include <tgmath.h>
#include <math.h>

#include "bessel.h"

#define I _Complex_I

#define NCHEB 50    /* order of the Chebyshev expansion */
#define NE 32        /* # of grid points for each dimension */
#define NX 1024       /* total # of grid points */
#define DIM 50        /* # of dimensions */
#define NSTEPS 600    /* # of propagation steps */
#define NDUMP 100     /* how often to report output */
#define dis 1     /* initial displacement */

/* Compute the stride length between indices to choose M indices
 almost uniformly from N options. N must be greater than or equal to M
 N >=M */

size_t uniform_stride(size_t N, size_t M)
{
    assert (N >= M);
    size_t stride = 1;
    size_t M1 = M-1;
    size_t N1 = N-1;
    while (stride * M1 < N1){
        stride++;
    }
    stride--;
    return stride;
}

static double * compute_Lp(size_t nx, double * nodes)
{
    double dx = nodes[1] - nodes[0];    /* Assumes constant dx. */
    for (size_t ii = 1; ii < nx-1; ii++){
        double dx2 = nodes[ii+1] - nodes[ii];
        if (fabs(dx2-dx) > 1e-15){
            fprintf(stderr, "Lp only works for uniform spacing\n");
            fprintf(stderr, "%3.15G %3.15G\n", dx, dx2);
            exit(1);
        }
            
    }
    
    double ub = nodes[nx-1] + dx; /* For periodic grid */
    double lb = nodes[0];

    double * Lp =calloc_double(nx  * nx);
    double dp = 2.0 * M_PI / (ub - lb);
    double * p = calloc_double(nx);
    for (size_t ii = 0; ii < nx; ii++){
        for (size_t jj = 0; jj < nx; jj++){
            Lp[ii*nx + jj] = 0.0;
        }
    }

    for (size_t ii = 0; ii < nx; ii++){
        p[ii] = dp * ii - dp * nx / 2.0;
    }

    for (size_t ll=0; ll<nx; ll++){
        for (size_t jj=0; jj<nx; jj++){
            for (size_t kk=0; kk<nx; kk++){
                double update =  creal(cexp(I*(nodes[jj]-nodes[ll])*p[kk])*pow(p[kk],2)*dx*dp/(2*M_PI));
                Lp[ll * nx + jj] = Lp[ll*nx + jj] - update;
            }
        }
    }
    free(p); p = NULL;

    return Lp;
}

struct GenericFunction * Lp_apply(const struct GenericFunction * f, void * args)
{
    double * Lp = args;
    assert (f->fc == LINELM);
    struct LinElemExp * ff = f->f;
    
    struct GenericFunction * g = generic_function_alloc(1, LINELM);

    size_t nx = ff->num_nodes;
    double * new_vals = calloc_double(nx);
    for (size_t jj = 0; jj < nx; jj++){
        new_vals[jj] = 0.0;
        for (size_t kk = 0; kk < nx; kk++){
            new_vals[jj] += Lp[kk*nx + jj] * ff->coeff[kk];
        }
    }
    g->f = lin_elem_exp_init(nx, ff->nodes, new_vals);
    free(new_vals); new_vals = NULL;
    return g;
}

/* initial state, given by the sum of a Gaussian state
and an additional term that aids sampling of high-dimensional wavefunctions
 when employing the the cross approximation (removed with function GSfix) */
double GS(const double *x)
{
  double out = exp(-0.5*pow(((x[0])-dis),2))/pow(M_PI,0.25);
  for(size_t j=1;j<DIM;j++){
    out *= exp(-0.5*pow((x[j]-dis),2))/pow(M_PI,0.25);
  }
  double out1 = exp(-0.5*pow((x[0]),2)/9.)/pow(M_PI,0.25);
  for(size_t j=1;j<DIM;j++){
    out1 *= exp(-0.5*pow((x[j]),2)/9.)/pow(M_PI,0.25);
  }
  out=1e28*(out+1e-5*out1);
  return out;
}

double GSfix(const double *x)
{
  double out = exp(-0.5*pow((x[0]),2)/9.)/pow(M_PI,0.25);
  for(size_t j=1;j<DIM;j++){
    out *= exp(-0.5*pow((x[j]),2)/9.)/pow(M_PI,0.25);
  }
  out=1e28*(1e-5*out);
  return out;
}

/* potential energy surface */
double Vpot(const double *x)
{
  double gamma = 0; /* DNA */
  double verticalscale=0.1;
  double out = verticalscale*(0.429*x[0]-1.126*pow(x[0],2)-0.143*pow(x[0],3)+0.563*pow(x[0],4));
  for(size_t jj=1;jj<DIM;jj++){
      out += verticalscale*(0.429*x[jj]-1.126*pow(x[jj],2)-0.143*pow(x[jj],3)+0.563*pow(x[jj],4));
      out += gamma*x[jj]*x[jj-1];
    }
  return out;
}

/* potential energy function for ft */
double Vpotft(const double *x, void * args)
{
    (void) args;
    double out = 0.0;
    out = Vpot(x);
    return out;
}


/* initial state gaussian for ft */
double psi0ft(const double *x, void * args)
{
    (void) args;
    double out = GS(x);
    return out;
}

/* term added to initial state gaussian for ft */
double psi0ftfix(const double *x, void * args)
{
    (void) args;
    double out = GSfix(x);
    return out;
}

void round_clean(struct FunctionTrain ** ft, double round_eps,
                 struct MultiApproxOpts * fopts)
{
    struct FunctionTrain * temp = function_train_round_maxrank_all(*ft, round_eps, fopts,5);
    function_train_free(*ft); *ft = NULL;
    *ft = function_train_copy(temp);
    function_train_free(temp); temp = NULL;
}

struct FunctionTrain * lap_prod_add(struct FunctionTrain * ft,
                                    double coeff, double DEm, double DEp,
                                    double m,
                                    struct FunctionTrain * V,
                                    double round_eps,
                                    struct Operator ** op,
                                    struct MultiApproxOpts * fopts)
{
  
    double * xtest = calloc_double(DIM);
    for (size_t l = 0; l<DIM; l++){
        xtest[l]=0.;
    }
  
    struct FunctionTrain * temp = exact_laplace_op(ft, op, fopts);
    function_train_scale(temp, -coeff/(DEm*m));
    round_clean(&temp, round_eps, fopts);
  
    struct FunctionTrain * temp2 = function_train_product(ft, V);
    function_train_scale(temp2, coeff * 2 / DEm);
    round_clean(&temp2, round_eps, fopts);

    struct FunctionTrain * out = function_train_copy(ft);
    function_train_scale(out, -coeff*DEp/DEm);
    c3axpy(1.0,temp, &out, 0.0, NULL);
    c3axpy(1.0,temp2, &out, 0.0, NULL);
    round_clean(&out, round_eps, fopts);

    function_train_free(temp);
    function_train_free(temp2);
    return out;
}

int clencheb(struct FunctionTrain **ft_psir,
           struct FunctionTrain **ft_psic,
           struct FunctionTrain *ft_cpolr[NCHEB],
           struct FunctionTrain *ft_cpoli[NCHEB],
           size_t N,
           double DEp, double DEm,
           double dtp, double dtm,
           struct FunctionTrain * ft_V, double m,
           struct Operator ** op,
           struct MultiApproxOpts * fopts )
{

    assert (N <= NCHEB);
    
    double round_eps=1e-8;
    size_t kk;
    double complex mik;
  
    /* For printing tests */
    double * xtest = calloc_double(DIM);
    for (size_t l = 0; l<DIM; l++){
        xtest[l]=0.;
    }
  
    for (size_t ii = 0; ii < 3; ii++){
        function_train_free(ft_cpolr[ii]);
        function_train_free(ft_cpoli[ii]);
        ft_cpolr[ii] = function_train_copy(*ft_psir);
        function_train_scale(ft_cpolr[ii], 0.);
        round_clean(&(ft_cpolr[ii]), round_eps, fopts);
        ft_cpoli[ii]=function_train_copy(*ft_psic);
        function_train_scale(ft_cpoli[ii], 0.);
        round_clean(&(ft_cpoli[ii]), round_eps, fopts);
      
    }
  
    for (size_t l=0; l<N; l++){
      
        ft_cpolr[2] = function_train_copy(ft_cpolr[1]);
        ft_cpoli[2] = function_train_copy(ft_cpoli[1]);
      
        ft_cpolr[1] = function_train_copy(ft_cpolr[0]);
        ft_cpoli[1] = function_train_copy(ft_cpoli[0]);
      
        kk = N-1-l;
        mik=pow(-I,kk);
      
        struct FunctionTrain * ft_temp1 = function_train_copy(*ft_psir);
        function_train_scale(ft_temp1, creal(mik));
        struct FunctionTrain * ft_temp2=function_train_copy(*ft_psic);
        function_train_scale(ft_temp2, -cimag(mik));
        c3axpy(1.0, ft_temp1, &ft_temp2, 0.0, NULL);
        function_train_scale(ft_temp2, bessj(kk, dtm));
        round_clean(&ft_temp2, round_eps, fopts);
        ft_cpolr[0] = function_train_copy(ft_temp2);
        c3axpy(-1.0, ft_cpolr[2], &(ft_cpolr[0]), 0.0, NULL);
        round_clean(&(ft_cpolr[0]), round_eps, fopts);

        struct FunctionTrain * ft_temp3 = function_train_copy(*ft_psir);
        function_train_scale(ft_temp3, cimag(mik));
        struct FunctionTrain * ft_temp4 = function_train_copy(*ft_psic);
        function_train_scale(ft_temp4, creal(mik));
        c3axpy(1.0, ft_temp3, &ft_temp4, 0.0, NULL);
        function_train_scale(ft_temp4, bessj(kk, dtm));
        round_clean(&ft_temp4, round_eps, fopts);
        ft_cpoli[0] = function_train_copy(ft_temp4);
        c3axpy(-1.0, ft_cpoli[2], &(ft_cpoli[0]), 0.0, NULL);
        round_clean(&(ft_cpoli[0]), round_eps, fopts);
      
        struct FunctionTrain * ft_temp5 = lap_prod_add(ft_cpolr[1], 2.0, DEm, DEp, m, ft_V, round_eps, op, fopts);
        struct FunctionTrain * ft_temp6 = lap_prod_add(ft_cpoli[1], 2.0, DEm, DEp, m, ft_V, round_eps, op, fopts);
      
        c3axpy(1.0, ft_temp5, &(ft_cpolr[0]), 0.0, NULL);
        round_clean(&(ft_cpolr[0]), round_eps, fopts);
        c3axpy(1.0, ft_temp6, &(ft_cpoli[0]), 0.0, NULL);
        round_clean(&(ft_cpoli[0]), round_eps, fopts);
      
        function_train_free(ft_temp1); ft_temp1 = NULL;
        function_train_free(ft_temp2); ft_temp2 = NULL;
        function_train_free(ft_temp3); ft_temp3 = NULL;
        function_train_free(ft_temp4); ft_temp4 = NULL;
        function_train_free(ft_temp5); ft_temp5 = NULL;
        function_train_free(ft_temp6); ft_temp6 = NULL;
      
        }
  
      c3axpy(-1.0, ft_cpolr[2], &(ft_cpolr[0]), 0.0, NULL);
      round_clean(&(ft_cpolr[0]), round_eps, fopts);
      function_train_free(*ft_psir); *ft_psir = NULL;
      *ft_psir = function_train_copy(ft_cpolr[0]);

  
      c3axpy(-1.0, ft_cpoli[2], &(ft_cpoli[0]), 0.0, NULL);
      round_clean(&(ft_cpoli[0]), round_eps, fopts);
      function_train_free(*ft_psic); *ft_psic = NULL;
      *ft_psic = function_train_copy(ft_cpoli[0]);

      struct FunctionTrain * ft_temp1 = function_train_copy(*ft_psir);
      function_train_scale(ft_temp1,cos(dtp));
      struct FunctionTrain * ft_temp2 = function_train_copy(*ft_psic);
      function_train_scale(ft_temp2,sin(dtp));
      c3axpy(1.0, ft_temp2, &ft_temp1, 0.0, NULL);

      struct FunctionTrain * ft_temp3 = function_train_copy(*ft_psic);
      function_train_scale(*ft_psic, cos(dtp));
      struct FunctionTrain * ft_temp4 = function_train_copy(*ft_psir);
      function_train_scale(*ft_psir,-sin(dtp));
      c3axpy(1.0, *ft_psir, ft_psic, 0.0, NULL);
      round_clean(ft_psic, round_eps, fopts);

      function_train_free(*ft_psir); *ft_psir = NULL;
      *ft_psir = function_train_round(ft_temp1, round_eps, fopts);
      
      function_train_free(ft_temp1); ft_temp1 = NULL;
      function_train_free(ft_temp2); ft_temp2 = NULL;
      function_train_free(ft_temp3); ft_temp3 = NULL;
      function_train_free(ft_temp4); ft_temp4 = NULL;

      for (size_t kk = 1; kk < N; kk++){
        if (kk % 10 == 0){
            printf("clenchebft kk = %zu\n", kk);
            printf("ft_psir_ranks ranks"); iprint_sz(DIM+1, function_train_get_ranks(*ft_psir));
            printf("ft_psic_ranks ranks"); iprint_sz(DIM+1, function_train_get_ranks(*ft_psic));
        }
    }
  
    return(0);
}

int main( int argc, char *argv[])
{
    (void) (argc);
    (void) (argv);

  double dt=0.01;              /* time increment */
    double m=1;                 /* mass */
    
    size_t dim = DIM;
    double lb = -5.;
    double ub = -lb;
    double time = 0;;
  
    /* beginning Chebyshe propagation */

    /* obtain psit at times dt*k by repeatedly applying the Nth-order Chebyshev expansion of exp(-I*H*dt) */

    /* grid parameters */
    double dx=(ub-lb)/NE;
    double dp=2*M_PI/(ub-lb);
    double xc[NE],pc[NE];
    double psier,psiei,psiem;
    double xo[DIM],po[DIM];

    size_t jy,jx;

    double Vmin= -1*DIM;
    double Vmax= 34.9455*DIM;
    double Emin=Vmin;
    double Emax=pow(M_PI,2)/(2*1*pow(dx,2))*dim+Vmax;
    double DEp=Emax+Emin;
    double DEm=Emax-Emin;
 
    double dtm=DEm*dt/2;
    double dtp=DEp*dt/2;
    double xir=0; /* real part of suvival ampliture */
    double xii=0; /* imaginary part of suvival ampliture */
    double normr=0; /* real part of norm */
    double normi=0; /* imaginary part of norm */
    double norm=0; /* total norm */

    /* grids */
    for (jy = 0; jy<NE; jy++){
        xc[jy] = lb + dx*jy;                /* coarse grid of coordinates */
        pc[jy] = dp*jy-dp*NE/2;                /* grid of momenta */
    }
    
    struct LinElemExpAopts * opts = lin_elem_exp_aopts_alloc(NE, xc);
    struct OneApproxOpts * qmopts = one_approx_opts_alloc(LINELM, opts);

    /* CROSS */
    struct C3Approx * c3a = c3approx_create(CROSS,dim);
    int verbose = 0;
    size_t init_rank = 3;
    double ** start = malloc_dd(dim);
    size_t stride = uniform_stride(NE, init_rank);
    struct Operator ** op = malloc(dim * sizeof(struct Operator *));
    assert (op != NULL);
    double * opLp[dim];
    for (size_t ii = 0; ii < dim; ii++){
        c3approx_set_approx_opts_dim(c3a,ii,qmopts);
        start[ii] = calloc_double(init_rank);
        for (size_t jj = 0; jj < init_rank; jj++){
            start[ii][jj] = xc[stride*jj];
        }

        opLp[ii] = compute_Lp(NE, xc);
        op[ii] = malloc(sizeof(struct Operator));
        op[ii]->f = Lp_apply;
        op[ii]->opts = opLp[ii];

    }
    c3approx_init_cross(c3a,init_rank,verbose,start);
    c3approx_set_adapt_kickrank(c3a,1);
    c3approx_set_cross_maxiter(c3a, 5);
    c3approx_set_cross_tol(c3a,1e-15);
    c3approx_set_round_tol(c3a,1e-15);
    c3approx_set_adapt_maxrank_all(c3a,10);
    free_dd(dim,start);
    struct MultiApproxOpts * fopts = c3approx_get_approx_args(c3a);

    /* FTs for the Vpot and psi0 built with previously defined functions  */
    size_t counter = 0;

    struct Fwrap * W_Vpotft = fwrap_create(dim,"general");
    fwrap_set_f(W_Vpotft,Vpotft,&counter);
    struct FunctionTrain * ft_V = c3approx_do_cross(c3a,W_Vpotft,1);
    
    struct Fwrap * W_psi0ft = fwrap_create(dim,"general");
    fwrap_set_f(W_psi0ft,psi0ft,&counter);
    struct FunctionTrain * ft_psir = c3approx_do_cross(c3a,W_psi0ft,1);
    struct FunctionTrain * ft_psic = function_train_constant(0, fopts);
    
    struct Fwrap * W_psi0ftfix = fwrap_create(dim,"general");
    fwrap_set_f(W_psi0ftfix,psi0ftfix,&counter);
    struct FunctionTrain * ft_psirfix = c3approx_do_cross(c3a,W_psi0ftfix,1);
    c3axpy(-1.0, ft_psirfix, &ft_psir, 0.0, NULL);
  
    double normctt=0; /* norm constant */
    normctt=function_train_inner(ft_psir,ft_psir);
    normctt=1/pow(normctt,0.5);
    function_train_scale(ft_psir, normctt);
  
    struct FunctionTrain * ft_psion = function_train_copy(ft_psir);


    struct FunctionTrain * ft_cpolr[NCHEB];
    struct FunctionTrain * ft_cpoli[NCHEB];
    for (size_t ii = 0; ii < NCHEB; ii++){
        ft_cpolr[ii] = NULL;
        ft_cpoli[ii] = NULL; 
    }

    printf("ft_psir ranks after cross"); iprint_sz(DIM+1, function_train_get_ranks(ft_psir));
    printf("ft_psir_ranks after cross"); iprint_sz(DIM+1, function_train_get_ranks(ft_psic));
    printf("ft_V_ranks after cross"); iprint_sz(DIM+1, function_train_get_ranks(ft_V));
    
    char append[20];
    char filename[256],fna[256],fnaa[256],fnaaa[256];
    size_t ncount=0;
    double * xx = calloc_double(dim);
    double v_ft_psir,v_ft_psic,v_ft_psim,ph;
    (void)v_ft_psir;
    (void)v_ft_psic;

    /* Propagation loop */

    strcpy(fna, "timings");
    FILE * fnn = fopen(fna,"w");
    
    strcpy(fnaa, "xi");
    FILE * fnx = fopen(fnaa,"w");
  
    strcpy(fnaaa, "norm");
    FILE * fny = fopen(fnaaa,"w");

    for(int k=0;k<=NSTEPS;k++){
        printf("k = %d/%u\n", k, NSTEPS);
        /* propagate */
        time=0;
        clock_t tic = clock();
        if(k>0){
            clencheb(&ft_psir,&ft_psic,ft_cpolr,ft_cpoli,NCHEB,DEp,DEm,dtp,dtm,ft_V,m,op,fopts);
        }
        xir=function_train_inner(ft_psir,ft_psion);
        xii=function_train_inner(ft_psic,ft_psion);
        fprintf(fnx,"%3.5G %3.5G %3.5G\n", dt*k, xir,xii);
        fflush(fnx);
      
        normr=function_train_inner(ft_psir,ft_psir);
        normi=function_train_inner(ft_psic,ft_psic);
        norm=(normr+normi);
        fprintf(fny,"%3.5G %3.5G\n", dt*k, norm);
        fflush(fny);

        clock_t toc = clock();
        time += (double)(toc - tic) / CLOCKS_PER_SEC;
        fprintf(fnn,"%d  %3.5G\n", k, time);
        fflush(fnn);

        if (k % NDUMP == 0){
            printf("k = %d\n", k);
            printf("ft_psir_ranks ranks"); iprint_sz(DIM+1, function_train_get_ranks(ft_psir));
            printf("ft_psir_ranks ranks"); iprint_sz(DIM+1, function_train_get_ranks(ft_psic));
            ncount += 1;
            strcpy(filename, "wave.");
            sprintf(append,"%zu",ncount);
            strcat(filename, append);
            printf("filename = %s\n", filename);
            FILE * fn = fopen(filename,"w");
            if(fn == NULL) {
                perror("fopen");
            }
            assert (fn != NULL);
            for (size_t j = 0; j<NX; j++){
                jy= j/NE;
                jx= j-jy*NE;

                for (size_t l = 0; l<DIM; l++){
                    xx[l]=0;
                }
                xx[0]=xc[jx];
                xx[1]=xc[jy];

                v_ft_psir=function_train_eval(ft_psir, xx);
                v_ft_psic=function_train_eval(ft_psic, xx);
                v_ft_psim=pow(v_ft_psir,2)+pow(v_ft_psic,2);

                fprintf(fn,"%3.5G %3.5G %3.5G\n",
                        xc[jx],xc[jy],v_ft_psim);

                if ((j+1) % NE == 0){
                    fprintf(fn,"\n");
                }
            }
            fclose(fn);
        }
    }

    /* cleanup memory */
    for (size_t ii = 0; ii < dim; ii++){
        free(op[ii]); op[ii] = NULL;
        free(opLp[ii]); opLp[ii] = NULL;
    }
    free(op); op = NULL;
    
    for (size_t ii = 0; ii < NCHEB; ii++){
        function_train_free(ft_cpolr[ii]); ft_cpolr[ii] = NULL;
        function_train_free(ft_cpoli[ii]); ft_cpoli[ii] = NULL;
    }
    fwrap_destroy(W_Vpotft); W_Vpotft = NULL; 
    fwrap_destroy(W_psi0ft); W_psi0ft = NULL;
    function_train_free(ft_psir); ft_psir = NULL;
    function_train_free(ft_psic); ft_psic = NULL;
    function_train_free(ft_V); ft_V = NULL;
    one_approx_opts_free_deep(&qmopts); qmopts = NULL;
    c3approx_destroy(c3a); c3a = NULL;

    return 0;
}
