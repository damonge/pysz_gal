import numpy as np
from ctypes import * 
import sys, os
from classy import Class

# sould match with that defined in FORTRAN codes?
MAXSIZE = 4096

# checking platform
LIBDIR = os.path.dirname(__file__)
# LIBDIR = (os.path.expanduser('~'))+'/Dropbox/py_mymodules/pysz'
sys.path.append(LIBDIR)


class tsz_gal_cl:
    def __init__(self):
        # print 'Class for tSZ Cl'
        # self.ptilde = np.loadtxt(LIBDIR+'/aux_files/ptilde.txt')
        self.fort_lib_cl = cdll.LoadLibrary(LIBDIR+"/source/calc_cl")

        self.fort_lib_cl.calc_cl_.argtypes = [
                                    POINTER(c_double), #h0
                                    POINTER(c_double), #obh2
                                    POINTER(c_double), #och2
                                    POINTER(c_double), #mnu
                                    POINTER(c_double), #bias
                                    POINTER(c_double), #Mcut
                                    POINTER(c_double), #M1
                                    POINTER(c_double), #kappa
                                    POINTER(c_double), #sigma_Ncen
                                    POINTER(c_double), #alp_Nsat
                                    POINTER(c_double), #rmax
                                    POINTER(c_double), #rgs
                                    POINTER(c_int64), #pk_nk
                                    POINTER(c_int64), #pk_nz
                                    np.ctypeslib.ndpointer(dtype=np.double), #karr
                                    np.ctypeslib.ndpointer(dtype=np.double), #pkarr
                                    np.ctypeslib.ndpointer(dtype=np.double), #dndz
                                    POINTER(c_int64), #nz_dndz
                                    POINTER(c_double), #z1
                                    POINTER(c_double), #z2
                                    POINTER(c_double), #z1_ng
                                    POINTER(c_double), #z2_ng
                                    POINTER(c_int64), #nl
                                    np.ctypeslib.ndpointer(dtype=np.double), #ell
                                    np.ctypeslib.ndpointer(dtype=np.double), #gg
                                    np.ctypeslib.ndpointer(dtype=np.double), #gy
                                    np.ctypeslib.ndpointer(dtype=np.double), #tll
                                    POINTER(c_double), #ng(z1<z<z2)
                                    POINTER(c_int64), #flag_nu
                                    POINTER(c_int64), #flag_tll
                                    POINTER(c_int64), #nm
                                    POINTER(c_int64) #nz
                                    ]
        self.fort_lib_cl.calc_cl_.restype = c_void_p

        # Calcualtion setup
        self.kmin = 1e-3
        self.kmax = 5.
        self.zmax = 4. # should be consistent with fortran code
        self.nk_pk = 200
        self.nz_pk = 51

        # Class
        self.cosmo = Class()

    def get_tsz_cl(self,ell_arr,params,dndz,z1,z2,z1_ng,z2_ng,nm,nz):
        self.zmin = z1
        self.zmax = z2
        obh2 = params['obh2']
        och2 = params['och2']
        As = params['As']
        ns = params['ns']
        mnu = params['mnu']
        mass_bias = params['mass_bias']
        Mcut = params['Mcut']
        M1 = params['M1']
        kappa = params['kappa']
        sigma_Ncen = params['sigma_Ncen']
        alp_Nsat = params['alp_Nsat']
        rmax = params['rmax']
        rgs = params['rgs']
        flag_nu_logic = params['flag_nu']
        flag_tll_logic = params['flag_tll']
        if type(flag_nu_logic) != bool:
            print 'flag_nu must be boolean.'
            sys.exit()
        if flag_nu_logic:
            flag_nu = 1
        else:
            flag_nu = 0
        if type(flag_tll_logic) != bool:
            print 'flag_tll must be boolean.'
            sys.exit()
        if flag_tll_logic:
            flag_tll = 1
        else:
            flag_tll = 0
        
        if 'theta' in params.keys():
            theta = params['theta']
            pars = {'output':'mPk','100*theta_s':theta,
                    'omega_b':obh2,'omega_cdm':och2,
                    'A_s':As,'n_s':ns,\
                    'N_ur':0.00641,'N_ncdm':1,'m_ncdm':mnu/3.,\
                    'T_ncdm':0.71611,\
                    'P_k_max_h/Mpc': self.kmax,'z_max_pk':self.zmax,\
                    'deg_ncdm':3.}
            self.cosmo.set(pars)
            self.cosmo.compute()
            h0 = self.cosmo.h()
        elif 'h0' in params.keys():
            h0 = params['h0']
            pars = {'output':'mPk','h':h0,
                    'omega_b':obh2,'omega_cdm':och2,
                    'A_s':As,'n_s':ns,\
                    'N_ur':0.00641,'N_ncdm':1,'m_ncdm':mnu/3.,\
                    'T_ncdm':0.71611,\
                    'P_k_max_h/Mpc': self.kmax,'z_max_pk':self.zmax,\
                    'deg_ncdm':3.}
            self.cosmo.set(pars)
            self.cosmo.compute()

        # get matter power spectra
        kh_arr = np.logspace(np.log10(self.kmin),np.log10(self.kmax),self.nk_pk)
        kh = np.zeros((self.nz_pk,self.nk_pk))
        pk = np.zeros((self.nz_pk,self.nk_pk))
        pk_zarr = np.linspace(self.zmin,self.zmax,self.nz_pk)
        for i in range(self.nz_pk):
            kh[i,:] = kh_arr
            if flag_nu == 0:
                pk[i,:] = np.array([self.cosmo.pk(k*h0,pk_zarr[i])*h0**3 for k in kh_arr])
            elif flag_nu == 1:
                pk[i,:] = np.array([self.cosmo.pk_cb(k*h0,pk_zarr[i])*h0**3 for k in kh_arr])

        # params
        h0_in = byref(c_double(h0))
        obh2_in = byref(c_double(obh2))
        och2_in = byref(c_double(och2))
        mnu_in = byref(c_double(mnu))
        mass_bias_in = byref(c_double(mass_bias))
        Mcut_in = byref(c_double(Mcut))
        M1_in = byref(c_double(M1))
        kappa_in = byref(c_double(kappa))
        sigma_Ncen_in = byref(c_double(sigma_Ncen))
        alp_Nsat_in = byref(c_double(alp_Nsat))
        rmax_in = byref(c_double(rmax))
        rgs_in = byref(c_double(rgs))
        flag_nu_in = byref(c_int64(flag_nu))
        flag_tll_in = byref(c_int64(flag_tll))

        # dNdz
        nz_dndz = byref(c_int64(len(dndz)))

        # integration setting
        z1_in = byref(c_double(self.zmin)) 
        z2_in = byref(c_double(self.zmax)) 
        
        # outputs
        nl = len(ell_arr)
        cl_gg = np.zeros((2,nl))
        cl_gy = np.zeros((2,nl))
        tll = np.zeros((nl*2,nl*2))
        ng = c_double(0.0)
        nl = c_int64(nl)
   
        self.fort_lib_cl.calc_cl_(
                h0_in, obh2_in, och2_in, mnu_in,\
                mass_bias_in, \
                Mcut_in, M1_in, kappa_in, sigma_Ncen_in, alp_Nsat_in,\
                rmax_in, rgs_in,\
                byref(c_int64(self.nk_pk)), byref(c_int64(self.nz_pk)),\
                np.array(kh),np.array(pk),\
                np.array(dndz),nz_dndz,\
                z1_in, z2_in,\
                byref(c_double(z1_ng)),byref(c_double(z2_ng)),\
                nl,np.array(ell_arr),\
                cl_gg,cl_gy,tll,ng,\
                flag_nu_in,flag_tll_in,\
                c_int64(nm), c_int64(nz)
                )

        self.cosmo.struct_cleanup()
        return cl_gg, cl_gy, tll, ng.value
