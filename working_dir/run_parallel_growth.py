import numpy as np
from chemdiff import *
from chemdiff.parallel import *
import os
import subprocess
from mpi4py import MPI

###### OPTIONAL #########
# if you don't have this module just comment out the `start` and `end` lines at beginning and end
from timeit import default_timer

start = default_timer()

################ INITIALIZE MPI ################
comm,nproc,rank = mpi_initialize()

############ SETUP AND CHEM PARAMS #############
f_pebout = 'pebble_composition.out'
with open(f_pebout,'w') as f:
	f.write('time    species    pebble_col_abundance\n')

grow_pebbles = True
# put your network file here
chm = '/home/ericvc/astrochem/networks/'
chm += 'data_chemnet_Deuterium_carbon_oxygen_frac.chm'


############ INITIALIZE THE COLUMN #############
r = 30*au
ti = 0
tf = 1.e4 # yrs

# temp profile from Krijt+2018
tmid = 130*(r/au)**(-1/2)
# surface density profile from Krijt+2018
Mdisk = msun*0.05
p = 1
rc = 100*au
sigc = (2-p)*Mdisk/(2*np.pi*rc*rc)
sig = sigc*(r/rc)**(-p) * np.exp(-(r/rc)**(2-p))

# model paramters
alpha = 1e-3
nzs = 50
dt = 100 # yrs
chemtime = 500 # yrs
touts = []
touts += [(i+1.)*pow(10,3) for i in range(9)]
touts += [(i+1.)*pow(10,4) for i in range(9)]
touts += [(i+2)/2*pow(10,5) for i in range(18)]
touts += [1.e6]
nts = len(touts)
nchems = int(tf/chemtime)
ndiffs = int(chemtime/dt)

# create the column
col = Column(r,tmid,alpha,nzs)
col.set_diff_params(dt)
h = col.h

############# INITIALIZE THE CELLS #############
o1618 = 500
o1617 = 2600

### Bosman 2018
init_abuns = {
	'H2' : 0.5,
	'He' : 9.75e-2,
	'NH3': 1.45e-6,
	'H2O': 1.18e-4,
	'CO' : 6.00e-5,
	'N2' : 2.00e-5,
	'CH4': 2.00e-6,
	'CH3OH' : 1.00e-6,
	'H2S' : 1.91e-8,
	'CO2' : 5.00e-5,
	'HCN' : 3.50e-7,
	'grain' : 2.2e-12
}
### Lyons and Young 2005
# init_abuns = {
# 	'H2' : 0.5,
# 	'He' : 0.16,
# 	'CO' : 2e-4,
# 	'H2O' : 2e-4,
# 	'grain' : 2.2e-12
# }
# ### isotopes
# init_abuns['H2-18-O'] = init_abuns['H2O']/o1618
# init_abuns['H2-17-O'] = init_abuns['H2O']/o1617
# init_abuns['C-18-O'] = init_abuns['CO']/o1618
# init_abuns['C-17-O'] = init_abuns['CO']/o1617

# Physical params
chi = 50
cosmic = 1.3e-17 # s^-1
grain_size = 0.1 # micron
dg0 = 0.01
opacity = 1000 # cm2 g-1 total opacity dust+gas
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
	tau += rho*100*opacity*col.dz*dg0

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
                              time,touts,chemtime=chemtime,network=chm)

	# rank 0 will do the pebble growth, diffusion, and update the cells
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
					peb_comp_norm = peb_comp[spec]/col.cells[0].NH
					f.write(f'{time}    {spec}    {peb_comp_norm:.10e}\n')

		goon = True
		# send the 'go on' signal and the updated column to other ranks
		for i in range(1,nproc):
			comm.send(True,dest=i,tag=15*i)
			comm.send(col,dest=i,tag=155*i)
	else:
		goon = comm.recv(source=0,tag=15*rank)
		col = comm.recv(source=0,tag=155*rank)

	while not goon:
		sleep(1)
		print('WAITING :: RANK ',rank)


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
