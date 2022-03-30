import numpy as np
from chemdiff import *
from chemdiff.parallel import *
import chemdiff.chemdiff_io as cdio
import os
import subprocess
from mpi4py import MPI

###### OPTIONAL #########
# if you don't have this module just comment out the `start` and `end` lines at beginning and end
from timeit import default_timer

start = default_timer()

################ INITIALIZE MPI ################
comm,nproc,rank = mpi_initialize()

############### READIN FROM FILE ################
model_inputs, phys_inputs, input_abundances = cdio.read_infile('cdinput.in')

############ SETUP AND CHEM PARAMS #############
f_pebout = model_inputs['pebfile']
with open(f_pebout,'w') as f:
	f.write('time    species    pebble_col_abundance\n')

grow_pebbles = True
# put your network file here
chm = model_inputs['chmfile']


############ INITIALIZE THE COLUMN #############
r = float(phys_inputs['r'])
if phys_inputs['r_units'] == 'au':
	r*=au
elif phys_inputs['r_units'] != 'cm':
	raise NameError("r_units must be au or cm, default is au")
ti = float(model_inputs['ti'])
tf = float(model_inputs['tf'])

# temp profile from Krijt+2018
tmid = 130*(r/au)**(-1/2)
# surface density profile from Krijt+2018
Mdisk = msun*0.05
p = 1
rc = 100*au
sigc = (2-p)*Mdisk/(2*np.pi*rc*rc)
sig = sigc*(r/rc)**(-p) * np.exp(-(r/rc)**(2-p))

# model paramters
alpha = float(phys_inputs['alpha'])
nzs = 50
dt = float(model_inputs['diff_dt'])
chemtime = float(model_inputs['chem_dt'])
touts = model_inputs['touts']
nts = len(touts)
nchems = int(tf/chemtime)

# create the column
col = Column(r,tmid,alpha,nzs)
h = col.h

############# INITIALIZE THE CELLS #############
o1618 = 500
o1617 = 2600

init_abuns = dict(input_abundances)

# Physical params
chi = float(phys_inputs['chi'])
cosmic = float(phys_inputs['cosmic'])
grain_size = float(phys_inputs['grain_size'])
dg0 = float(phys_inputs['dg0'])
opacity = float(phys_inputs['opacity'])
rho0 = sig/np.sqrt(2.*np.pi)/h
xray = 0
zq = 3
tatm = 2*tmid

if rank == 0:
	if not os.path.exists('r00'):
		os.system('mkdir r00')

# set column densities to zero
NCO = 0.
NH2 = 0.
NH = 0.
tau = 0.
for j in reversed(range(nzs)):
	dirr=f'r00/z{j:0>2}'
	if rank == 0:
		if not os.path.exists(dirr):
			os.system('mkdir '+dirr)
	z = col.dz*(j+0.5) # cm
	rho = rho0*np.exp(-z*z/2./h/h)
	# temp profile from krijt+2018
	temp = tmid+(tatm-tmid)*(np.sin(np.pi*z/2./zq/h))**(4)
	if z >= zq*h:
		temp = tatm
	# column densities
	nh = 2*rho/mbar
	nco = 0
	if 'CO' in init_abuns:
		nco = init_abuns['CO']*nh
	nh2 = 0
	if 'H2' in init_abuns:
		nh2 = init_abuns['H2']*nh
	NCO += nco*col.dz
	NH2 += nh2*col.dz
	NH += nh*col.dz
	# optical depth
	tau += rho*opacity*col.dz*dg0

	### CREATE THE CELL
	col.cells[j] = Cell(r,z,chi=chi,cosmic=cosmic,grain_size=grain_size,dust_gas_ratio=dg0,
		av=tau/3.02,rho=rho,Tgas=temp,Tdust=temp,xray=xray,NCO=NCO,NH2=NH2,NH=NH,
		abundances = dict(init_abuns))

################ MAIN LOOP #####################
peb_comp = {}
cwd = os.getcwd()
for t in range(nchems):
	time = (t+1)*chemtime

	# do chemistry
	do_parallel_chemistry(col,comm,nproc,rank,
                              time,touts,verbose=True,chemtime=chemtime,network=chm)


##################### TIMER ####################
if rank == 0:
	end = default_timer()
	time_to_run = end-start
	hrs = time_to_run//3600
	mins = (time_to_run-(hrs*3600))//60
	secs = time_to_run-(hrs*3600+mins*60)
	with open('time_to_run.out','a') as f:
		f.write(f'time to run : {time_to_run:.1f} sec\n')
		f.write(f'            = {int(hrs)}h {int(mins)}m {secs:.1f}s\n')

print('done!')
