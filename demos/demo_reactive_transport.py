"""
Reactive transport modeling using xGems
"""
#%%
from __future__ import division,print_function
import sys,os
myxgemspath = '../build/lib' #change this to path where xgems is compiled. If directory structure is same as in repo dont change this
sys.path.append(myxgemspath)
from xgems import ChemicalEngine
import numpy as np
from scipy.special import erfc
import matplotlib.pylab as plt
from copy import deepcopy as copy
import time
#set matplotlib font size and type
font = {'family' : 'serif',
         'size'   : 12}
plt.rc('font',**font)
#%%
class VarDef(object):
    """
    a class to define variable used for convenience 
    """
    _args_name=['type','default','doc']
    def __init__(self,*args):
        for name,val in zip(self._args_name,args):
            setattr(self,name,val)
#%%
class xGems(object):
    """
    xGems object for coupling with reactive transport
    """
    _vars={'gems_init_file':VarDef('str','','Gems file for initialization ends with .lst'),
           'gems_ic_files':VarDef('list',[],'Gems dbr files list for setting initial conditions'),
           'ic_cell_labels':VarDef('scalar',0,'an array specify cell label corresponding to initial state given by dbr file'),
           'gems_bc_files':VarDef('list',[],'Gems dbr files for setting boundary conditions'),
           'T':VarDef('scalar',0,'temprature for each cell'),
           'pH':VarDef('scalar',0,'pH'),
           'pe':VarDef('scalar',0,'pe'),
           'P':VarDef('scalar',0,'pressure for each cell'),
           'nx':VarDef('param',0,'number of cells'),
            '_poros':VarDef('scalar',1,'internal variable for porosity'),
            '_phaseVolFrac':VarDef('dict',1,'internal for  phase volume fractions'),
            '_phaseConc':VarDef('dict',1,'internal for  phase concentration'),
            'pH':VarDef('scalar',0,'pH'),
            'pE':VarDef('scalar',0,'pe'),
            }
    
    def __init__(self,domain,inputs={}):
        """
        initalizes gems reaction module
        
        Input
        -----
         inputs: dict
             dictionary containing inputs for initializing xGems classs
             inputs['gems_init_file']: input file obtained from GEMS3K to initialize xGems
             inputs['gems_ic_files']: list containing files obtained from GEMS3K for initial conditions
             inputs['gems_bc_files']:  list containing files obtained from GEMS3K for boundary conditions
             inputs['ic_cell_labels']: array to specify label for cell corresponding to index of files in initial conditions
        """
        self.nx = domain.nx  
        self._read_inputs(inputs)  
        self.gem = ChemicalEngine(self.gems_init_file) 
        self.nelements = self.gem.numElements()
        self.element_names = []
        for i in range(self.nelements):
            self.element_names.append(self.gem.elementName(i))
        self.nphases= self.gem.numPhases()
        self.phase_names =[]
        for i in range(self.nphases):
            self.phase_names.append(self.gem.phaseName(i))
        
    def _read_inputs(self,inputs):
        """
        a convenience method to read inputs
        """
        skip_items=['nx']
        for k,v in self._vars.iteritems():
            if k not in skip_items:
                if v.type=='scalar':
                    setattr(self,k,inputs.get(k,v.default)*np.array([1]*self.nx))
                else:
                    setattr(self,k,inputs.get(k,v.default))
            
    def advance(self,c_dict):
        """
        advances a timestep in calculation
        
        Input
        -----
        c_dict: dict
            dictionary of concentrations of elements obtained from  transport step
        """
        self.c_mob = copy(c_dict)
        for i in range(self.nx):
            T,P = self.T[i],self.P[i]
            b=[]
            for name in self.element_names:
                b.append(c_dict[name][i]+self.c_immob[name][i])
            self.gem.equilibrate(T,P,b)
            T,P,b = self.gem.temperature(),self.gem.pressure(),self.gem.elementAmounts()
            c_mob  = self.gem.elementAmountsInPhase(self.phase_names.index('aq_gen'))
            for j,name in enumerate(self.element_names):
                self.c_mob[name][i]=c_mob[j]
                self.c_immob[name][i]=b[j]-c_mob[j]
            pVol = self.gem.phaseVolumes()
            pAmt  = self.gem.phaseAmounts()
            V = self.gem.systemVolume()
            for j,name in enumerate(self.phase_names):
                self._phVfrac[name][i] = pVol[j]/V
                self._phConc[name][i] = pAmt[j]             
            self._poros[i]=self._phVfrac['aq_gen'][i]
            self.pH[i] = self.gem.pH()
            self.pE[i]=self.gem.pe()
        return 
 
    def get_ic(self):
        """
        provides element concentrations for initializing transport solver 
        
        Returns
        -------
        c_dict: dict
            dictionary of concentrations of elements  
        """
        out=[]
        for fname in self.gems_ic_files:
            T,P,b=self.read_lst_file(fname)
            self.gem.equilibrate(T, P, b)  
            T = self.gem.temperature() 
            P = self.gem.pressure()  
            b = self.gem.elementAmounts()  
            V = self.gem.systemVolume()
            c_mob  = self.gem.elementAmountsInPhase(self.phase_names.index('aq_gen'))
            ph,pe =self.gem.pH(),self.gem.pe() 
            pVol =self.gem.phaseVolumes()
            pAmount = self.gem.phaseAmounts()
            out.append([T,P,V,c_mob,b-c_mob,ph,pe,pVol,pAmount])
        c_dict={}
        self.c_immob ={}
        self.c_mob = {}
        for name in self.element_names:
            c_dict[name]=[]
            self.c_mob[name]=[]  
            self.c_immob[name]=[]
        self._phConc ={}
        self._phVfrac={}
        self._poros = []
        for name in self.phase_names:
            self._phConc[name]=[]
            self._phVfrac[name]=[]            
        for i in xrange(self.nx):
            idx = self.ic_cell_labels[i]
            T,P,V,c_mob,c_immob,ph,pe,phaseVol,phaseAmt = out[idx]
            self.T[i]=T
            self.P[i]=P
            for j,name in enumerate(self.element_names):
                self.c_immob[name].append(c_immob[j])
                self.c_mob[name].append(c_mob[j])
                c_dict[name].append(c_mob[j])  
            for j,name in enumerate(self.phase_names):
                self._phConc[name].append(phaseAmt[j])             
                self._phVfrac[name].append(phaseVol[j]/V)
                if name == 'aq_gen':
                    self._poros.append(phaseVol[j]/V)
            self.pH[i] = ph
            self.pE[i]=pe
        return c_dict
    
    def get_bc(self):
        """
        provides element concentration for 
        """
        out=[]
        idx = 0
        for fname in self.gems_bc_files:
            T,P,b=self.read_lst_file(fname)
            self.gem.equilibrate(T, P, b)
            T = self.gem.temperature() 
            P = self.gem.pressure() 
            b = self.gem.elementAmounts() 
            c_mob  = self.gem.elementAmountsInPhase(self.phase_names.index('aq_gen'))
            out.append({})
            for i,name in enumerate(self.element_names):
                out[idx][name]=c_mob[i]
            idx+=1
        return out

    def get_cdict(self):
        """
        returns element concentration for 
        """
        return self.c_mob

    @property
    def porosity(self):
        """
        porosity
        """
        return self._poros 
   
    @property
    def phase_volume_frac(self):
        """
        phase volume fractions
        """
        return self._phVfrac
    
    @property
    def phase_conc(self):
        """
        phase concentrations
        """
        return self._phConc
    
    @staticmethod
    def read_lst_file(fname):
        gem = ChemicalEngine(fname)
        return gem.temperature(),gem.pressure(),gem.elementAmounts()
    
#%%
class Domain1D(object):
    """
    Initializes 1D domain (cell-centered)
    """
    def __init__(self,L,dx):
        self.dx = dx
        self.L = L
        self.x = np.linspace(self.dx/2, self.L-self.dx/2, round(self.L/self.dx))
        self.nx = len(self.x)
        

#%%
class Ade1D(object):
    """
    Explicit cell-centered finite difference solver for advection-diffusion equation. Uses upwinding scheme for advection and
    central difference for diffusion.
    ::math
       partial_{t}C + \partial_{x}J = 0
       J = -D0*a*phi^n \partial_{x} C
    """
    _vars={'L':VarDef('param',0,'length of domain'),
           'dx':VarDef('param',0,'grid spacing (only uniform grid supported)'),
           'D0':VarDef('param',0,'diffusion coefficient in pore space'),
           'dt':VarDef('param',0,'timestep'),
           't':VarDef('param',0,'time elapsed'),
           'c':VarDef('scalar',0,'concentration'),
           'oldc':VarDef('scalar',0,'concentration from previous timestep'),
           'u':VarDef('scalar',0,'velocity'),
           'x':VarDef('scalar',0,'node co-ordinates'), 
           'phi':VarDef('scalar',1,'porosity'),
           'bc':VarDef('dict',{'left':['c',0],'right':['j',0]},'definition of boundary conditions'),
           'archie_params':VarDef('dict',{'a':1,'n':1},'archie relationship parameters'),
          }

    def __init__(self,domain,inputs={}):
        """
        intializes ADE1D
        """
        self._read_domain(domain)
        self._read_inputs(inputs)
        self.dt =  0.1*self.dx**2/np.max(self.Dp)
        
    def _read_domain(self,domain):
        """
        reads domain
        """
        self.dx,self.L,self.x,self.nx = domain.dx,domain.L,domain.x,domain.nx
        
    def _read_inputs(self,inputs):
        """
        reads input
        """
        skip_items=['L','dx','x','nx']
        for k,v in self._vars.iteritems():
            if k not in skip_items:
                if v.type=='scalar':
                    setattr(self,k,inputs.get(k,v.default)*np.array([1]*self.nx))
                else:
                    setattr(self,k,inputs.get(k,v.default))
    
    def advance(self):
        """
        advances one time step through explicit time stepping
        """
        self.oldc = copy(self.c)
        c =self.c
        Dp,dx,u =self.Dp,self.dx,self.u
        dt = 0.5*self.dx**2/np.max(Dp)
        self.dt = dt
        self.t += dt
        De_r,De_l = np.roll(Dp,-1), np.roll(Dp, 1)
        De_l,De_r = (De_l+Dp)/2,(De_r+Dp)/2
        cr,cl = np.roll(c,-1),np.roll(c, 1)
        self.c = c - u * (c-cl)*(dt/(dx)) * (u<=0) - u * (cr-c)*(dt/dx) * (u>0) \
        + dt/(dx**2)*(De_l*cl+De_r*cr-(De_l+De_r)*c)  
        self.apply_bc()
        
    def apply_bc(self):
        """
        applies boundary conditions
        """
        Dp,c,dx=self.Dp,self.c,self.dx    
        def setC(cl,cr):
            return (cl+cr)/2
        #left bc
        if self.bc['left'][0]=='c':
            self.c[0]=setC(self.bc['left'][1],c[1])
        if self.bc['left'][0]=='j':
            j=self.bc['left'][1]
            cl = c[1]+(2*dx*j/Dp[0])
            self.c[0]=setC(cl,c[1])
        #right bc
        if self.bc['right'][0]=='c':
            c[-1]=setC(self.bc['right'][1],self.c[-2])
        if self.bc['right'][0]=='j':
            j=self.bc['right'][1]
            cr = c[-2]-(2*dx*j/Dp[0])
            self.c[-1]=setC(c[-2],cr)
            
    @property
    def Dp(self):
        """
        pore diffusion coefficient
        """
        a,n,D0,phi=self.archie_params['a'],self.archie_params['n'],self.D0,self.phi
        return D0*a*phi**(n-1)

    def run(self,tf,iout=100,verbose=True):
        """
        runs solver till time equals to tf
        """
        i=0
        t0=time.time()
        while self.t <tf:
            i+=1
            if verbose==True:
                if (i%iout)==0: print ('time:',self.t)
            self.advance()
            if verbose==True:
                print ('time:',self.t)
                print ('+++Time taken for simulation (in s):',time.time()-t0)
                
#%%
class MultiComponentTransport(object):
    """
    class for multicomponent transport assuming same diffusion coefficient for all species
    """
    def __init__(self,domain,inputs,component_list):
        """
        intializes 
        """
        self.dx,self.L,self.x,self.nx = domain.dx,domain.L,domain.x,domain.nx
        self.component_list = component_list
        self.components={}
        for name in component_list:
            local_inputs = copy(inputs)
            if name in local_inputs:
                local_inputs.update(local_inputs[name])
            self.components[name] = Ade1D(domain,local_inputs)
            
    def get_cdict(self):
        """
        get concentration of all species in a dictionary 
        """
        out= {}
        for k,v in self.components.iteritems():
            out[k] =copy(v.c)
        return out
        
    def advance(self):
        """
        advances one time step
        """
        for k,v in self.components.iteritems():
            v.advance()
            
    def run(self,tf,iout=100,verbose=True):
        """
        runs solver till t equals to tf
        """
        i=0
        t0=time.time()
        while self.t <tf:
            i+=1
            if verbose==True:
                if (i%iout)==0: print ('time:',self.t)
            self.advance()
            if verbose==True:
                print ('time:',self.t)
                print ('+++Time taken for simulation (in s):',time.time()-t0)
                
    def update_c(self,c_dict):
        """
        updates concentration field
        """
        for k,v in self.components.iteritems():
            setattr(v,'c',copy(c_dict[k]))
    
    @property
    def Dp(self):
        """
        pore diffusion coefficient
        """
        return self.components[self.component_list[0]].Dp
    
    @property
    def t(self):
        """
        current time
        """
        return self.components[self.component_list[0]].t

    @property
    def phi(self):
        """
        porosity
        """
        return self.components[self.component_list[0]].phi

    @property
    def dt(self):
        """
        timestep
        """
        return self.components[self.component_list[0]].dt
#%%
def test_diffusion():
    """
    tests diffusion implementation
    """
    def analytical(cb,D,x,tf):
        """
        Analytical solution for 1D advection diffusion equation
        """
        c = cb*erfc((x)/(4*D*tf)**0.5) 
        return c
    inputs = {}
    cb=0.01
    inputs['bc']={'left':['c',cb],'right':['c',0]}
    inputs['D0']=1
    domain = Domain1D(50,2)
    model = Ade1D(domain,inputs)
    model.run(30,verbose=False)
    c= analytical(cb,model.Dp,model.x,model.t)
    err =np.sqrt(np.mean((model.c-c)**2)/np.mean(c**2))
    assert err < 5e-2
    print ("relative error:%s"%err )
    if True:
        plt.plot(model.x,model.c, label='simulated')
        plt.plot(model.x,c,'o',label='analytical')
        plt.legend()
        plt.xlabel('distance')
        plt.ylabel('concentration')
        plt.show()
    return model
#%%    
def test_advection_diffusion():
    """
    tests Advection diffusion implementation
    """
    def analytical(cb,x,ts,D, u):
        c = (cb/2)*(erfc((x - u*ts)/(4*D*ts)**0.5)+
             erfc((x+ u * ts)/(4*D *ts)**0.5)*np.exp(u*x/D))
        return c
    inputs = {}
    cb=0.01
    inputs['bc']={'left':['c',cb],'right':['c',0]}
    inputs['D0']=2.1e-8
    inputs['u']=1e-5
    domain = Domain1D(0.5,0.5/800.)
    model = Ade1D(domain,inputs)
    model.run(model.dt*2000,verbose=False)
    c= analytical(cb,model.x,model.t,model.Dp,model.u)
    err =np.sqrt(np.mean((model.c-c)**2)/np.mean(c**2))
    assert err < 5e-2
    print ("relative error:%s"%err )
    if True:
        plt.plot(model.x,model.c, label='simulated')
        plt.plot(model.x,c,'--',label='analytical')
        plt.legend()
        plt.xlabel('distance')
        plt.ylabel('concentration')
        plt.show()
    return model
#%%
def test_multicomponent_transport():
    """
    tests multicomponent reactive transport implementation
    """
    def analytical(cb,D,x,tf):
        """
        Analytical solution for 1D advection diffusion equation
        """
        c = cb*erfc((x)/(4*D*tf)**0.5) 
        return c
    inputs = {}
    cb_Na=0.01
    cb_Cl = 0.02
    inputs['Na'] = {}
    inputs['Na']['bc']={'left':['c',cb_Na],'right':['c',0]}
    inputs['Cl']={}
    inputs['Cl']['bc']={'left':['c',cb_Cl],'right':['c',0]}
    inputs['D0']=1
    domain = Domain1D(50,2)
    model = MultiComponentTransport(domain,inputs,['Na','Cl'])
    model.run(30,verbose=False)
    c_Na= analytical(cb_Na,model.Dp,model.x,model.t)
    c_Cl= analytical(cb_Cl,model.Dp,model.x,model.t)
    if True:
        plt.plot(model.x,model.components['Na'].c, label='Na simulated')
        plt.plot(model.x,c_Na,'o',label='Na analytical')
        plt.plot(model.x,model.components['Cl'].c, label='Cl simulated')
        plt.plot(model.x,c_Cl,'o',label='Cl analytical')
        plt.legend()
        plt.xlabel('distance')
        plt.ylabel('concentration')
        plt.show()
    return model

#%%        
def test_xgems_reactive_trnasport():
    """
    tests gems reactive transport for benchmark problem
    ref: Shao H et al. Appl. Geochem. 2009 Jul 1, 24(7):1287-300.
    """
    #reactive transport example
    domain=Domain1D(0.5,0.5/80) #set up simulation domain
    inputs={}
    inputs['phi']=0.32
    inputs['u']= 9.375e-6
    inputs['D0'] = 9.375e-6 * 0.0067
    inputs['gems_init_file']=os.path.join('resources','CalciteIC-dat.lst')
    inputs['gems_ic_files']=[os.path.join('resources','CalciteIC-dat.lst')]
    inputs['gems_bc_files']=[os.path.join('resources','CalciteBC-dat.lst')]
    inputs['ic_cell_labels'] = 0
    rxn= xGems(domain,inputs) #initalize xgems
    #get initial and boundary conditions from gems
    ic = rxn.get_ic()
    bc = rxn.get_bc()
    #setup correct_boundary and initial conditions
    for name in rxn.element_names:
        inputs[name]={}
        inputs[name]['c'] = ic[name]
        inputs[name]['bc']= {'left':['c',bc[0][name]],'right':['j',0]}
    trans= MultiComponentTransport(domain,inputs,rxn.element_names) #initialize transport
#    run model we use sequential non-interative approach for coupling 
    i=0
    while trans.t < 21000:
        i+=1
        trans.advance() #run transport 
        cdict = trans.get_cdict() #get conc
        rxn.advance(cdict) #supply conc and run gems
        cdict = rxn.get_cdict() #get new equilibriated concentrations
        trans.update_c(cdict) #update concentration in transpor module
        if i%100==0:print ("Time: %s"%trans.t)
#   post processing of results
    fig,ax1= plt.subplots()
    cdict = trans.get_cdict()
    x = trans.x
    outlist=['Ca','Mg','Cl']
    lineformat = ['-k','-b','-r','-g','-c']
    i=0
    for name in outlist:
        ax1.plot(x,cdict[name],lineformat[i],label = name)
        i+=1
    ax1.set_ylabel('ionic concentration (M)')
    ax1.set_xlabel('Distance (m)')
    ax2 = ax1.twinx()
    outlist=['Calcite','Dolomite-dis']
    for name in outlist:
        ax2.plot(x,rxn.phase_conc[name],lineformat[i],label = name)
        i+=1
    ax2.set_ylabel('mineral concentration (M)')
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc=1,fontsize=10)   
    ax2.set_ylim([0,4e-4])
    ax1.ticklabel_format(style='sci',axis='y',scilimits=(0,0),useMathText=True)
    ax2.ticklabel_format(style='sci',axis='y',scilimits=(0,0),useMathText=True)
    plt.show()        
    return trans,rxn
    
if __name__ == '__main__':
    trans,rxn=test_xgems_reactive_trnasport()
