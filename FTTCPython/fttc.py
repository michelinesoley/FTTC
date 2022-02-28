import numpy as np
from scipy.special import jv
import c3py

def GS(r,param=None):
    global dis
    nevals,dim=r.shape
    out=np.zeros((nevals,))
    for ii in range(dim):
      if ii==0:
        out=np.exp(-0.5*(r[:,ii]-dis)**2)/np.pi**0.25
      else:
        #out=out*np.exp(-0.5*(r[:,ii])**2)/np.pi**0.25 # not displaced
        out=out*np.exp(-0.5*(r[:,ii]-dis)**2)/np.pi**0.25 # displaced
    out1=np.zeros((nevals,))
    for ii in range(dim):
      if ii==0:
        out1=np.exp(-0.5*(r[:,ii])**2/9.)/np.pi**0.25
      else:
        out1=out1*np.exp(-0.5*(r[:,ii])**2/9.)/np.pi**0.25
    out=1e7*(out+1e-5*out1)
    return out

def GSfix(r,param=None):
    nevals,dim=r.shape
    out1=np.zeros((nevals,))
    for ii in range(dim):
      if ii==0:
        out1=np.exp(-0.5*(r[:,ii])**2/9.)/np.pi**0.25
      else:
        out1=out1*np.exp(-0.5*(r[:,ii])**2/9.)/np.pi**0.25
    out1=1e7*(1e-5*out1)
    return out1

def V(r,param=None):
    gamma=0.
    verticalscale=0.1
    nevals,dim=r.shape
    out=np.zeros((nevals,))
    for ii in range(dim):
        if ii == 0: # Double well potential
          out=out+verticalscale*(0.429*r[:,ii]-1.126*r[:,ii]**2-0.143*r[:,ii]**3+0.563*r[:,ii]**4)
        else:
          out=out+verticalscale*(0.429*r[:,ii]-1.126*r[:,ii]**2-0.143*r[:,ii]**3+0.563*r[:,ii]**4)
          out=out+verticalscale*gamma*r[:,ii]*r[:,ii-1]
#      out=out+0.5*(r[:,ii])**2 # Harmonic oscillator potential
    return out
    
def zeros(r,param=None):
    nevals,dim=r.shape
    out=np.zeros((nevals,))
    return out
    
def lap_prod_add(ft,coeff,dem,dep,m,ft_V,dim,lb,ub,nparam,basis,op):
    global round_eps,rmax
    temp=ft.laplace_op(op)
    temp.scale(-coeff/(dem*m))
    temp.round(eps=round_eps,maxrank_all=rmax)
    temp2=ft*ft_V
    temp2.scale(coeff*2/dem)
    temp2.round(eps=round_eps,maxrank_all=rmax)
    out = ft.copy()
    out=out.scale(-coeff*dep/dem)
    out=temp+out
    out.round(eps=round_eps,maxrank_all=rmax)
    out=temp2+out
    out.round(eps=round_eps,maxrank_all=rmax)
    return out
    
def clencheb(ft_psir,ft_psic,ncheb,dep,dem,dtp,dtm,ft_V,m,dim,lb,ub,nparam,basis,op): # Clenshaw Chebyshev propagation

    global round_eps,rmax
    
    ft_cpolr=[None]*3
    ft_cpoli=[None]*3
    
    for ii in range(3):
    
      ft_cpolr[ii]=build_ft(zeros,None,dim,lb,ub,nparam,basis=basis)
      ft_cpolr[ii].round(eps=round_eps,maxrank_all=rmax)
      ft_cpoli[ii]=build_ft(zeros,None,dim,lb,ub,nparam,basis=basis)
      ft_cpoli[ii].round(eps=round_eps,maxrank_all=rmax)

    for l in range(ncheb):
    
      ft_cpolr[2]=ft_cpolr[1].copy()
      ft_cpoli[2]=ft_cpoli[1].copy()
      ft_cpolr[1]=ft_cpolr[0].copy()
      ft_cpoli[1]=ft_cpoli[0].copy()
      
      kk=ncheb-1-l
      mik=(-complex(0,1))**kk
      
      ft_temp1=ft_psir.copy()
      ft_temp1.scale(np.real(mik))
      ft_temp2=ft_psic.copy()
      ft_temp2.scale(-np.imag(mik))
      ft_temp2=ft_temp1+ft_temp2
      ft_temp2.round(eps=round_eps,maxrank_all=rmax)
      ft_temp2.scale(jv(kk,dtm))
      ft_temp2.round(eps=round_eps,maxrank_all=rmax)
      ft_cpolr[0]=ft_temp2.copy()
      ft_cpolr[0]=ft_cpolr[0]-ft_cpolr[2]
      ft_cpolr[0].round(eps=round_eps,maxrank_all=rmax)

      ft_temp3=ft_psir.copy()
      ft_temp3.scale(np.imag(mik))
      ft_temp4=ft_psic.copy()
      ft_temp4.scale(np.real(mik))
      ft_temp4=ft_temp3+ft_temp4
      ft_temp4.round(eps=round_eps,maxrank_all=rmax)
      ft_temp4.scale(jv(kk,dtm))
      ft_temp4.round(eps=round_eps,maxrank_all=rmax)
      ft_cpoli[0]=ft_temp4.copy()
      ft_cpoli[0]=ft_cpoli[0]-ft_cpoli[2]
      ft_cpoli[0].round(eps=round_eps,maxrank_all=rmax)
      
      ft_temp5=lap_prod_add(ft_cpolr[1],2.0,dem,dep,m,ft_V,dim,lb,ub,nparam,basis,op)
      ft_temp6=lap_prod_add(ft_cpoli[1],2.0,dem,dep,m,ft_V,dim,lb,ub,nparam,basis,op)
      
      ft_cpolr[0]=ft_temp5+ft_cpolr[0]
      ft_cpolr[0].round(eps=round_eps,maxrank_all=rmax)
      ft_cpoli[0]=ft_temp6+ft_cpoli[0]
      ft_cpoli[0].round(eps=round_eps,maxrank_all=rmax)
          
    ft_cpolr[0]=ft_cpolr[0]-ft_cpolr[2]
    ft_cpolr[0].round(eps=round_eps,maxrank_all=rmax)
    ft_psir=ft_cpolr[0].copy()
    
    ft_cpoli[0]=ft_cpoli[0]-ft_cpoli[2]
    ft_cpoli[0].round(eps=round_eps,maxrank_all=rmax)
    ft_psic=ft_cpoli[0].copy()
    
    ft_temp1=ft_psir.copy()
    ft_temp1.scale(np.cos(dtp))
    ft_temp2=ft_psic.copy()
    ft_temp2.scale(np.sin(dtp))
    ft_temp1=ft_temp2+ft_temp1
    ft_temp1.round(eps=round_eps,maxrank_all=rmax)
    
    ft_temp3=ft_psic.copy()
    ft_psic.scale(np.cos(dtp))
    ft_temp4=ft_psir.copy()
    ft_psir.scale(-np.sin(dtp))
    ft_psic=ft_psir+ft_psic
    ft_psic.round(eps=round_eps,maxrank_all=rmax)
 
    ft_psir=ft_temp1.copy()
    ft_psir.round(eps=round_eps,maxrank_all=rmax)
      
    return ft_psir,ft_psic

def build_ft(func,args,dim,lb,ub,nparam,basis):

    ft = c3py.FunctionTrain(dim)
    for ii in range(dim):
        if basis == "piecewise":
            ft.set_dim_opts(ii,"piecewise",lb[ii],ub[ii],nparam=nparam,coeff_check=1,tol=1e-1,nregions=2)
        elif basis == "polynomial":
            ft.set_dim_opts(ii,"legendre",lb[ii],ub[ii],nparam=nparam,maxnum=40,coeff_check=1,tol=1e-3)
        elif basis == "linelm":
            ft.set_dim_opts(ii,"linelm",lb[ii],ub[ii],nparam=nparam,)

    verbose=0
    init_rank=3
    adapt=0

    maxrank=10
    kickrank=1
    crosstol=1e-15
    roundtol=1e-15
    maxiter=5

    ft.build_approximation(func,args,init_rank,verbose,
                           adapt,maxrank=maxrank,
                           round_tol=roundtol,
                           cross_tol=crosstol,
                           kickrank=kickrank,
                           maxiter=maxiter)
    return ft

# Dynamics algorithm

def run():

    global round_eps,rmax
    global dis

    ncheb=50 # order of Chebyshev expansion
    ne=32 # number of grid points per dimension
    nx=ne**2 # total number of grid points
    dim=47 # number of dimensions
    nsteps=600 # number of propagation steps
    ndump=100 # how often to report output
    dis=1 # initial displacement
    dt=0.01 # time increment
    m=1 # mass
    xmin=-5 # lower bound
    xmax=5 # upper bound
    dx=(xmax-xmin)/ne # positiongrid division
    dp=2*np.pi/(xmax-xmin) # momentum grid division
    vmin=-1*dim#(-0.0849)*dim # minimum potential # double well potential
#    vmin=0 # Harmonic oscillator potential
    vmax=34.9455*dim#33.9455*dim # maximum potential # double well potential
#    vmax=0.5*(xmax**2)*dim
    emin=vmin # minimum energy
    emax=np.pi**2/(2*dx**2)*dim+vmax # maximum energy
    dep=emax+emin
    dem=emax-emin
    dtm=dem*dt/2
    dtp=dep*dt/2
    round_eps=1e-8
    rmax=5
    
    # define grids
    xc=np.zeros(ne)
    pc=np.zeros(ne)
    for jy in range(ne):
      xc[jy]=xmin+dx*jy
      pc[jy]=dp*jy-dp*ne/2
    op = c3py.c3.build_lp_operator(dim,xc) # (ne,xc)
    
    # define function trains
    basis="linelm"
    nparam=ne # number of linear elements
    # support of univariate functions
    lb=[xc[0]]*dim #[xmin]*dim
    ub=[xc[-1]]*dim #[xmax]*dim

    # build function trains
    emptyline=np.asarray([[]])
    ft_V=build_ft(V,None,dim,lb,ub,nparam,basis=basis)
    ft_V.round(eps=round_eps,maxrank_all=rmax)
    potfile=open("potential.npy","w")
    xcur=np.zeros(dim)
    for jj in range(nx):
      jy=int(jj/ne)
      jx=int(jj-jy*ne)
      for ll in range(dim):
        xcur[ll]=0
      xcur[0]=xc[jx]
      xcur[1]=xc[jy]
      potentialvalue=ft_V.eval(xcur)
      potfileline=np.asarray([[xc[jx],xc[jy],potentialvalue]])
      np.savetxt(potfile,potfileline)
      if (jj+1)%ne==0:
        np.savetxt(potfile,emptyline)
    potfile.close()
    
    ft_psir=build_ft(GS,None,dim,lb,ub,nparam,basis=basis)
    ft_psir.round(eps=round_eps,maxrank_all=rmax)
    ft_psirfix=build_ft(GSfix,None,dim,lb,ub,nparam,basis=basis)
    ft_psirfix.round(eps=round_eps,maxrank_all=rmax)
    ft_psir=ft_psir-ft_psirfix
    normcur=ft_psir.inner(ft_psir)
    ft_psir.round(eps=round_eps,maxrank_all=rmax)
    ft_psir.scale(1/np.sqrt(normcur))
    ft_psir.round(eps=round_eps,maxrank_all=rmax)
    ft_psic=build_ft(zeros,None,dim,lb,ub,nparam,basis=basis)
    ft_psic.round(eps=round_eps,maxrank_all=rmax)
    ft_psion=ft_psir.copy() # initial state
    ft_psion.round(eps=round_eps,maxrank_all=rmax)
    
    # initialize arrays
    xir=np.zeros(nsteps+1)
    xii=np.zeros(nsteps+1)
    norm=np.zeros(nsteps+1)
    print("Step No. = ",0)
    xir[0]=ft_psir.inner(ft_psion)
    xii[0]=ft_psic.inner(ft_psion)
    print("Autocorrelation Function = ",xir[0],xii[0])
    normr=ft_psir.inner(ft_psir)
    normi=ft_psic.inner(ft_psic)
    norm[0]=normr+normi
    print("Norm = ",norm[0])
    
    xifile=open("xi.npy","w")
    normfile=open("norm.npy","w")
    xifileline=np.asarray([[0.,xir[0],xii[0]]])
    np.savetxt(xifile,xifileline)
    xifile.flush()
    normfileline=np.asarray([[0.,norm[0]]])
    np.savetxt(normfile,normfileline)
    normfile.flush()
    wavefile=open("wave00.npy","w")
    xcur=np.zeros(dim)
    for jj in range(nx):
      jy=int(jj/ne)
      jx=int(jj-jy*ne)
      for ll in range(dim):
        xcur[ll]=0
      xcur[0]=xc[jx]
      xcur[1]=xc[jy]
      v_ft_psir=ft_psir.eval(xcur)
      v_ft_psic=ft_psic.eval(xcur)
      v_ft_psim=v_ft_psir**2+v_ft_psic**2
      wavefileline=np.asarray([[xc[jx],xc[jy],v_ft_psim]])
      np.savetxt(wavefile,wavefileline)
      if (jj+1)%ne==0:
        np.savetxt(wavefile,emptyline)
    wavefile.close()
    
    # Propagate
  
    dumpcount=0
    for ii in range(1,nsteps+1):
      print("Step No. = ",ii)
      ft_psir,ft_psic=clencheb(ft_psir,ft_psic,ncheb,dep,dem,dtp,dtm,ft_V,m,dim,lb,ub,nparam,basis,op)
      xir[ii]=ft_psir.inner(ft_psion)
      xii[ii]=ft_psic.inner(ft_psion)
      print("Autocorrelation Function = ",xir[ii],xii[ii])
      normr=ft_psir.inner(ft_psir)
      normi=ft_psic.inner(ft_psic)
      norm[ii]=normr+normi
      print("Norm = ",norm[ii])
      dumpcount=dumpcount+1
      xifileline=np.asarray([[ii*dt,xir[ii],xii[ii]]])
      np.savetxt(xifile,xifileline)
      xifile.flush()
      normfileline=np.asarray([[ii*dt,norm[ii]]])
      np.savetxt(normfile,normfileline)
      normfile.flush()
      if ii%ndump==0:
        wavefile=open("wave{0}.npy".format(dumpcount),"w")
        xcur=np.zeros(dim)
        for jj in range(nx):
          jy=int(jj/ne)
          jx=int(jj-jy*ne)
          for ll in range(dim):
            xcur[ll]=0
          xcur[0]=xc[jx]
          xcur[1]=xc[jy]
          v_ft_psir=ft_psir.eval(xcur)
          v_ft_psic=ft_psic.eval(xcur)
          v_ft_psim=v_ft_psir**2+v_ft_psic**2
          wavefileline=np.asarray([[xc[jx],xc[jy],v_ft_psim]])
          np.savetxt(wavefile,wavefileline)
          if (jj+1)%ne==0:
            np.savetxt(wavefile,emptyline)
        wavefile.close()
    xifile.close()
    normfile.close()

if __name__ == "__main__":

    run()
