import numpy as np
import matplotlib.pyplot as plt
import astrochem.wrapper as ac
import chemdiff as cd
from chemdiff.constants import *
from chemdiff.chemdiff_io import get_defaults

def temp_profile(r):
	return 130*(r/au)**(-1/2)

def surface_density_profile(r):
	Mdisk = msun*0.05
	p = 1
	rc = 100*au
	sigc = (2-p)*Mdisk/(2*np.pi*rc*rc)
	return sigc*(r/rc)**(-p) * np.exp(-(r/rc)**(2-p))

r = 30*au

tmid = temp_profile(r)
sig = surface_density_profile(r)

## model parameters
tf = 1.e6
model_defaults, phys_defaults, abun_defaults = get_defaults()
alpha = 1e-3
nzs = 50
dt = 100 # yrs
chemtime = 100 # yrs
touts = model_defaults['touts']
nts = len(touts)
nchems = int(tf/chemtime)
ndiffs = int(chemtime/dt)

## create the column
col = cd.Column(r,tmid,alpha,nzs)
col.set_diff_params(dt)
h = col.h

#### initialize cells
init_abuns = dict(abun_defaults)

