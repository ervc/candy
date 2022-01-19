import matplotlib.pyplot as plt
import numpy as np
from chemdiff import *
import os
import subprocess
from mpi4py import MPI
from timeit import default_timer
from time import sleep

start = default_timer()

def mpi_initialize():
    comm = MPI.COMM_WORLD
    nproc = comm.size
    rank = comm.rank
    return comm,nproc,rank

def do_parallel_chemistry(col,size,time,touts,chemtime=500,
    network='network.chm',abs_err=1e-20,rel_err=1e-10):
    '''
    Runs Astrochem in parallel

    PARAMETERS
    ----------
    col : Column object to run astrochem on
    size : int, number of processers to use (nproc)
    time : float, the current time, check if output should be saved
    touts : list, times to save outputs for
    astrochem_params: args passed to run astrochem
        chemtime, network, abs_err, rel_err
    '''
    # print('starting chem, rank',rank)
    my_cells = None
    nzs = col.ncells
    js_per_rank = nzs//size
    remain = nzs - (size*js_per_rank)
    if rank == 0:
        for i in range(1,size):
            comm.send(col.cells[i*js_per_rank+remain:(i+1)*js_per_rank+remain],dest=i,tag=i*11)
        my_cells = col.cells[0:js_per_rank+remain]
    else:
        my_cells = comm.recv(source=0,tag=rank*11)

    cwd = os.getcwd()
    wait = True
    if rank == 0:
        all_done = np.zeros(size,dtype='bool')
    for my_j,cell in enumerate(my_cells):
        # get the absolute j for each 
        if rank == 0:
            j = my_j
        else:
            j = rank*js_per_rank+remain+my_j
        # print(rank, j)
        dirr = f'{cwd}/r00/z{j:0>2}'
        cell.write_chem_inputs(chemtime,abs_err=abs_err,rel_err=rel_err,
            f_net=network,f_input=dirr+'/input.ini',f_source=dirr+'/source.mdl')

        subprocess.run(['astrochem','-q','input.ini'],cwd=dirr)
        # print(' ran ',rank,j)
        if float(time) in touts:
            subprocess.run(['cp','astrochem_output.h5',f'astrochem_output_t{time:0>7}.h5'],cwd=dirr)
            subprocess.run(['cp','source.mdl',f'source_t{time:0>7}.mdl'],cwd=dirr)

        # update cell abundances with chem abundances
        out_times,d = get_abundict(f'{dirr}/astrochem_output.h5','all')
        for spec in d:
            cell.abundances[spec] = d[spec][-1]
    # print('finished cells, rank',rank)

    # once all the cells have done send a signal to rank0
    # I'm sure there's a better way to do this, but I find the back and forth
    # send and receive to/from all the ranks is a good way to make sure
    # everyone is done with what they're supposed to do

    # In this first one, all ranks other than 0 send a signal they are done
    # if rank zero is first here, it will wait for all the signals
    if rank != 0:
        comm.send(True,dest=0,tag=rank*12)
    else:
        all_done[0] = True
        for i in range(1,size):
            all_done[i] = comm.recv(source=i,tag=i*12)

    # then here rank 0 says "Ok, everyone is good to go here"
    # other ranks will wait here for the go ahead from rank 0
    if rank == 0:
        if all(all_done):
            wait = False
            for i in range(1,size):
                comm.isend(False,dest=i,tag=i*13)
    else:
        wait = comm.recv(source=0,tag=13*rank)

    # by the time we get to here, no one should be waiting, but
    # this remains as a good check to make sure everyone is together
    while wait:
        print('waiting :: rank ',rank)
        sleep(1)
    # print('finished chem, rank',rank)


def grow_grains(col,peb_comp,time,grow_pebbles=True):
    '''
    Grows pebbles near the midplane with timescale from Birnsteil 2012

    PARAMETERS
    ----------
    col : Column object
    peb_comp : dict, current pebble composition dictionary
    time : float, current time
    grow_pebbles : bool, should pebbles be grown or not?

    OUTPUT
    ------
    dict, updated peb_comp dictionary
    '''
    # print('growing grains,rank', rank)
    nzs = col.ncells
    for j in range(nzs):
        cell = col.cells[j]
        # change grain abundances if z < h
        t_grow = 1./(cell.dust_gas_ratio*col.omega)
        if cell.z/col.h <= 1 and time > 0 and grow_pebbles:
            deps = -col.dt*sec/(t_grow + time*sec)
        else:
            deps = 0
        for spec in cell.abundances:
            if spec[:5] == 'grain':
                d_ice = deps*cell.abundances[spec] # X*dt/(tau+t) = X*f*dt
                cell.abundances[spec] += d_ice
                if spec not in peb_comp:
                    peb_comp[spec] = 0
                peb_comp[spec] -= d_ice*cell.nh*col.dz # col averaged abundance
    return peb_comp


def do_diffusion(col):
    '''
    Do the diffusion calculation for the column

    PARAMETERS
    ----------
    col : Column Object
    '''
    # print('doing diffusion, rank',rank)
    nzs = col.ncells
    col_abunds = col.get_abundance_array()
    newarray = {}
    for spec in col_abunds:
        newarray[spec] = np.zeros(nzs)
        for j in range(nzs):
            rhos_j = col.cells[j].rho
            sp = col_abunds[spec]
            fp = 0.
            if j < nzs-1:
                rhos_j1 = col.cells[j+1].rho
                fp = 0.5*(rhos_j+rhos_j1)*col.beta*(sp[j+1]-sp[j])/col.dz
            elif j == nzs-1:
                fp = 0.
            if j == 0:
                fm = -fp
            newarray[spec][j] = sp[j]+(fp-fm)/col.dz/rhos_j
            fm = fp
    
    # save new diffusion values
    for spec in newarray:
        for j in range(nzs):
            col.cells[j].dust_gas_ratio = grain_abun2dg(newarray['grain'][j])
            col.cells[j].abundances[spec] = newarray[spec][j]

    
def update_cells(col,opacity=1000):
    '''
    Update the column densities and optical depths

    PARAMETERS
    ----------
    col : Column Object
    '''
    print('updating cells,rank',rank)
    nzs = col.ncells

    NCO = 0.
    NH2 = 0.
    NH = 0.
    tau = 0.
    for j in reversed(range(nzs)):
        cell = col.cells[j]
        nh = cell.nh
        nco = 0
        if 'CO' in list(cell.abundances.keys()):
            nco = cell.abundances['CO'] * nh
        nh2 = 0
        if 'H2' in list(cell.abundances.keys()):
            nh2 = cell.abundances['H2'] * nh
        NCO += nco * col.dz
        NH2 += nh2 * col.dz
        NH += nh * col.dz # this should stay the same
        tau += cell.rho*opacity*col.dz*100*cell.dust_gas_ratio
        cell.NCO = NCO
        cell.NH2 = NH2
        cell.NH = NH
        cell.av = tau/3.02
    


###### Initialize MPI ######
comm,nproc,rank = mpi_initialize()
############################

######################################## SET UP ###############################################

f_pebout = 'pebble_composition.out'
with open(f_pebout,'w') as f:
    f.write('time    species    pebble_col_abundance\n')

grow_pebbles = True

chm = '/Users/ericvc/astrochem/networks/'
chm += 'umist12_x.chm'

####################################### INITIALIZE COLUMN ########################################

r = 100*au
tf = 1.e4 # yrs

# temp profile from Krijt 2018
tmid = 130*(r/au)**(-1/2) # K 
# surface density profile from Krijt 2018
Mdisk = msun*0.05
p = 1
rc = 100*au
sigc = (2-p)*Mdisk/(2*np.pi*rc*rc)
sig = sigc * (r/rc)**(-p) * np.exp(-(r/rc)**(2-p))

# model parameters
alpha = 1e-3
nzs = 50
dt = 100 # yrs
chemtime = 500 # yrs
touts = []
touts += [(i+1.)*pow(10,3) for i in range(9)]
touts += [(i+1.)*pow(10,4) for i in range(9)]
touts += [(i+2.)/2*pow(10,5) for i in range(18)]
touts += [1.e6] # set up log tout list
nts = len(touts)
nchems = int(tf/chemtime)
ndiffs = int(chemtime/dt)

#create column
col = Column(r,tmid,alpha,nzs)
col.set_diff_params(dt)
h = col.h

############################### INITIALIZE CELLS #############################

init_abuns = {
    'H2' : 0.5,
    'He' : 9.75e-2,
    'NH3': 1.45e-6,
    'H2O': 1.18e-4,
    'CO' : 6.00e-5,
    'N2' : 2.00e-5,
    'CH4': 2.00e-6,
    'CH3OH' : 1.00e-6,
    'H2S': 1.91e-8,
    'CO2': 5.00e-5,
    'HCN' : 3.5e-7,
    'grain' : 2.2e-12
}

# phys params
chi = 50
cosmic = 1.3e-17 # s-1
grain_size = 0.1 # micron
dg0 = 0.01
opacity = 1000 # cm2 g-1
rho0 = sig/np.sqrt(2.*np.pi)/h
xray = 0
zq = 3
tatm = 2*tmid

#set up cells in column
if not os.path.exists('r00'):
    os.system('mkdir r00')

# initialize column densities as zero
NCO = 0.
NH2 = 0.
NH = 0.
tau = 0.
av1 = 49
# note work from top down => reversed(range(nzs))
# Initialize the cells for all ranks
for j in reversed(range(nzs)):
    dirr = f'r00/z{j:0>2}'
    if not os.path.exists(dirr):
            os.system('mkdir '+dirr)
    z = col.dz*(j+0.5) # cm

    #hydrostatic equilibrium for density
    rho = rho0*np.exp(-z*z/2./h/h)
    # temp profile from krijt 2018
    temp=tmid+(tatm-tmid)*(np.sin(np.pi*z/2./zq/h))**(4)
    if z >= zq*h:
            temp = tatm

    # column densities
    nh = 2*rho/mbar
    nco = 0
    if 'CO' in init_abuns:
        nco = init_abuns['CO'] * nh
    nh2 = 0
    if 'H2' in init_abuns:
        nh2 = init_abuns['H2'] * nh
    NCO += nco * col.dz
    NH2 += nh2 * col.dz
    NH += nh * col.dz

    # optical depth
    tau += rho*opacity*col.dz*100*dg0
    if tau/3.02 <= 1.0:
        av1 = j
    col.cells[j] = Cell(r,z,chi=chi,cosmic=cosmic,grain_size=grain_size,dust_gas_ratio=dg0,
        av=tau/3.02,rho=rho,Tgas=temp,Tdust=temp,xray=xray,NCO=NCO,NH2=NH2,NH=NH,
        abundances=dict(init_abuns))

all_abunds = {}
zs = [c.z/col.h for c in col.cells]


# plot temp, rho, av
fig,axs = plt.subplots(1,3,figsize=(15,5))

ax = axs[0]
ax.plot(zs,[c.Tgas for c in col.cells])
ax.set(xlabel='Z [Scale Height]',ylabel='Temp [K]')

ax = axs[1]
ax.plot(zs,[c.rho for c in col.cells])
ax.set(xlabel='Z [Scale Height]',ylabel='Density [g cm$^{-3}$]',yscale='log')

ax = axs[2]
ax.plot(zs,[c.av for c in col.cells])
ax.set(xlabel='Z [Scale Height]',ylabel='Av',yscale='log')

plt.savefig('t0_phys_params.png')
plt.close()

################################## MAIN LOOP #####################################

peb_comp = {}
cwd = os.getcwd()
print(cwd)
for t in range(nchems):

    time = (t+1)*chemtime
    
    # do chem on each cell
    do_parallel_chemistry(col,nproc,time,touts,chemtime=chemtime,network=chm)
    
    # rank 0 will do the pebble growth, diffusion, and update the cells
    # this part is not too long anyway
    goon = False
    if rank == 0:
        for diff_loop in range(ndiffs):
            # change grain abundances
            peb_comp = grow_grains(col,peb_comp,time,grow_pebbles)
            do_diffusion(col)
        update_cells(col,opacity)

        if float(time) in touts:
            print(time)
            col_abunds = col.get_abundance_array()
            with open(f_pebout,'a') as f:
                for spec in peb_comp:
                    peb_comp_norm = peb_comp[spec]/col.cells[0].NH # record composition relative to column abundance
                    f.write(f'{time}    {spec}    {peb_comp_norm:.5e}\n')

        goon = True
        for i in range(1,nproc):
            comm.send(True,dest=i,tag=15*i)
    else:
        goon = comm.recv(source=0,tag=15*rank)

    while not goon:
        sleep(1)
        print('waiting :: rank ',rank)

        
if rank == 0:
    end = default_timer()

    time_to_run = end-start
    hrs = time_to_run//3600
    mins = (time_to_run-(hrs*3600))//60
    secs = time_to_run-(hrs*3600+mins*60)
    with open('time_to_run.out','a') as f:
        f.write(f'time to run : {time_to_run:.1f} sec\n')
        f.write(f'            = {int(hrs)}h {int(mins)}m {secs:.1f}s\n')

print('done')


