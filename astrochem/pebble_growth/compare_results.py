import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import chemdiff as cd

nzs = 50
ts = np.arange(1,11)
ts *= int(1e3)
nts = ts.size

specs = ['CO','H2O','C','O']
cmap = cm.get_cmap('magma',len(specs)+1)

for t,time in enumerate(ts):
	print(t)
	serials = {}
	para2s = {}
	para3s = {}
	zs = np.zeros(nzs)
	for j in range(nzs):
		z = j*(5/nzs)+(5/nzs/2)
		zs[j] = z

		ser_dir = f'serial/r00/z{j:0>2}'
		par_dir = f'parallel/para2/r00/z{j:0>2}'
		par3_dir = f'parallel/para3/r00/z{j:0>2}'

		# serial
		_,d = cd.get_abundict(ser_dir+f'/astrochem_output_t{time:0>7}.h5',specs)
		for spec in specs:
			if spec not in serials:
				serials[spec] = np.zeros(nzs)
			serials[spec][j] = d[spec][-1]

		# parallel 2
		_,d = cd.get_abundict(par_dir+f'/astrochem_output_t{time:0>7}.h5',specs)
		for spec in specs:
			if spec not in para2s:
				para2s[spec] = np.zeros(nzs)
			para2s[spec][j] = d[spec][-1]

		# parallel 3
		_,d = cd.get_abundict(par3_dir+f'/astrochem_output_t{time:0>7}.h5',specs)
		for spec in specs:
			if spec not in para3s:
				para3s[spec] = np.zeros(nzs)
			para3s[spec][j] = d[spec][-1]

	# print(serials['CO'])
	fig,ax = plt.subplots()
	for k,spec in enumerate(specs):
		ax.plot(zs,serials[spec],ls='-',c=cmap(k),label=spec)
		ax.plot(zs,para2s[spec],ls='--',c=cmap(k))
		ax.plot(zs,para3s[spec],ls=':',c=cmap(k))
	ax.set(xscale='linear',yscale='log',xlabel='Z/H',ylabel='Abundance')
	ax.set(title=f't = {time} yrs')
	ax.plot([],[],ls='-',c='k',label='Serial')
	ax.plot([],[],ls='--',c='k',label='2 Cores')
	ax.plot([],[],ls=':',c='k',label='3 Cores')
	ax.legend()

	plt.savefig(f'compare_abund_{t:0>2}.png')
	plt.close()