

import pylab as plt
import argparse
from sys import stdout
import numpy as np
from scipy import interpolate
import matplotlib.colors as colors
import seaborn as sb
from bisect import bisect_left,bisect_right

print('libs imported')

pars=argparse.ArgumentParser()
pars.add_argument('-cf','--configFile',type = str, help = "File containing configuration arguments and values")
pars.add_argument('-g','--gfile',type = str,help = "File containing genes")
pars.add_argument('-i','--infile',type = str,help = "Bedgraph file containing the signal")
pars.add_argument('-gt','--gfiletype',type = str,help = "Type of the file containing genes ['bed','scoretsv']; default = 'bed'")
pars.add_argument('-gi','--gindx',nargs = 5,type = str)
pars.add_argument('-rp','--replot',type = str,help = "File containing matrix to be replotted")

pars.add_argument('-oh','--hmfile',type = str,help = "Heatmap output filename")
pars.add_argument('-oa','--avgfile',type = str,help = "Average profile output filename")
pars.add_argument('-om','--matfile',type = str,help = "matrix output filename")

pars.add_argument('-r','--region',type = str, help = "Region to be plotted ['TSS','TSE','genebody']; default = 'TSS'")
pars.add_argument('-fl','--flank',type = str, help = "Length of flanking fragments to be plotted with the selected region; default = 3000")

pars.add_argument('-only','--chrmonly',nargs = '+',type = str, help = "The exact names of chromosomes to be  exclusively considered")
pars.add_argument('-co','--chrmomit',nargs = '+',type = str, help = "The exact names of chromosomes to be  exclusively considered")
pars.add_argument('-go','--gomit',nargs = '+',type = str,help = "Names of features to be ignored")
pars.add_argument('-nb','--nbest',type = str,help = "The number of features with the best score to be used")
pars.add_argument('-sr','--scorerange',nargs = 2,type = str,help = "The score range from which the features are selected")
pars.add_argument('-of','--ofirst',action = 'store_true', help = "Whether to use only the first feature with the same name [no value]; default = False")

pars.add_argument('-nt','--nticks',type = str)
pars.add_argument('-p','--plottype',type = str,help = "Type of plot to be generated ['avgprof','heatmap','both']; default = avgprof")
pars.add_argument('-s','--sort',action = 'store_true', help = "Whether to sort the matrix used for generating a heatmap [no value]; default = False")
pars.add_argument('-cm','--cmap',type = str,help = "Colormap used in the heatmap; default = 'Reds'")
pars.add_argument('-sm','--smooth',type = str,help = "Smoothing factor used when smoothing the average profile with a spline. Set to 'false' or 0 if you don't want to smooth; default = flank*1e-4")
pars.add_argument('-ht','--hmtitle',type = str,help = "Title of the heatmap")
pars.add_argument('-at','--avgtitle',type = str,help = "Title of the average profile")
pars.add_argument('-cb','--cbar',action = 'store_true',help = 'Whether to show a colorbar next to the heatmap [no value]; default = False')
pars.add_argument('-hn','--hnorm',type = str,help = "Type of norm to be used for the heatmap colorscale ['lin','log']; default = 'lin'")
args=pars.parse_args()


config = {		#default argument values are stored here
'region':'TSS',
'flank':'1000',

'infile':None,
'gfile':None,
'gfiletype':'bed',
'gindx':None,
'replot':None,

'hmfile':None,
'avgfile':None,
'matfile':None,

'plottype':'avgprof',
'cmap':'Reds',
'nticks':'1',
'sort':False,
'smooth':False,
'hmtitle':None,
'avgtitle':None,
'hnorm':'lin',
'cbar':False,

'chrmomit':None,
'chrmadd':None,
'gomit':None,
'chrmonly':None,

'nbest':None,
'scorerange':None,
'ofirst':False
}

print('args parsed')

def readConfig():
	if args.configFile:
		f = open(args.configFile,'r')

		for line in f:
			l = line.split()
			if len(l)==2:
				if l[0] in config:
					config[l[0]] = l[1]
		f.close()
	
	#override the config file values with the console values
	for arg in vars(args):
		if getattr(args, arg):
			config[arg] = getattr(args, arg)

	for arg in ['sort','ofirst','cbar']:
		if type(config[arg]) == str:
			config[arg] = config[arg].lower() == 'true'
	for arg in ['gomit','chrmomit','gindx','scorerange','chrmonly']:
		if type(config[arg]) == str:
			config[arg] = config[arg].split(',')
	for arg in ['flank','nticks','nbest']:
		if type(config[arg]) == str:
			config[arg] = int(config[arg])
	for arg in ['region','plottype','hnorm','gfiletype']:
		if type(config[arg]) == str:
			config[arg] = config[arg].lower()

def loadGenes():

	c = {}
	if config['chrmonly']:
		for x in config['chrmonly']:
			c[x] = []

	# Add a new file format here
	# read<Format>() should return a list of iterables with format:
	# (<chromosome_name>,<TSS>,<TSE>,<additional_info>...)
	if config['gfiletype'] == 'bed':
		genes = readBed()
	elif config['gfiletype'] == 'scoretsv':
		genes = readScore()

	for g in genes:
		if not config['chrmonly'] and g[0] not in c:
			c[g[0]] = []
		if g[0] in c:
			c[g[0]].append((g[1],g[2]))				

	keys = list(c.keys())
	for key in c:
		if region == 'genebody':
			c[key] = sorted(c[key],key = lambda g: min(g))
		else:
			c[key] = sorted(c[key],key = lambda g: g[TSE])
	print('loadGenes done')
	return c

def query(fname,genome):

	f = open(fname,'r')
	result = []
	line = f.readline().split()
	
	while line:
		chrm = line[0]
		if chrm in genome:

			
			signal = np.zeros(int(line[1]),dtype = 'int8')

			
			start = 0
			end = int(line[1])
			

			i = 1
			for gene in genome[chrm]:
				stdout.write(chrm+'\t'+str(i)+'\r')
				i+=1
				### Rozszerzenie szukanego framentu o fragmenty flankujące
				if region == 'genebody':
					gStart = min(gene) - flank
					gEnd = max(gene) + flank +1
				else:
					gStart = gene[TSE] - flank
					gEnd = gene[TSE] + flank +1
			
				if gStart >= end:
					while line and int(line[1]) < gStart:
						if line[0] != chrm:
							break
						x,y,z = int(line[1]),int(line[2]),int(line[3])
						line = f.readline().split()
					if gStart < y:
						signal = np.append(np.array([z]*(y-gStart),dtype = 'int8'),np.zeros(int(line[1])-y,dtype = 'int8'))
					else:
						signal = np.zeros(int(line[1])-gStart,dtype = 'int8')
					start = gStart
					end = int(line[1])

				elif gStart > start:
					signal = signal[gStart-start:]	
					start = gStart

				elif gStart < start and start == 0:
					signal = np.append(np.zeros(start - gStart,dtype = 'int8'),signal)
					start = gStart


				if gEnd+1 > end:
					signal = np.append(signal,np.zeros(gEnd-end,dtype = 'int8'))
					x,y,z = int(line[1]),int(line[2]),int(line[3])
					signal[x-start:y-start] = z
					end = y
					line = f.readline().split()
					while gEnd+1 > end:
						if not line or line[0] != chrm:
							end = gEnd
							break
						x,y,z = int(line[1]),int(line[2]),int(line[3])
						signal[x-start:y-start] = z
						end = y
						line = f.readline().split()
					if gEnd < x:
						signal = np.append(signal,np.zeros(x-gEnd,dtype = 'int8'))
					if gEnd < y:
						signal = np.append(signal,np.array([z]*(y-x),dtype = 'int8'))
			

				if gene[0] > gene[1]:
					res = np.flip(signal[:gEnd-start],axis = 0)
				else:
					res = signal[:gEnd-start]

				if region == "genebody":
					res = np.append(res[:flank],[normalize(res[flank:-flank],flank),res[-flank:]])
				
				result.append(res)
			
			print(chrm,'done',i-1)

			while line and line[0] == chrm:
				line = f.readline().split()				
		else:
			while line and line[0] == chrm:
				line = f.readline().split()	
	f.close()
	return np.array(result,dtype = 'int8')

def normalize(arr,size):

	spl = interpolate.InterpolatedUnivariateSpline(np.linspace(0,len(arr),len(arr)),arr,k=3)
	return spl(np.linspace(0,len(arr),size))

def plot(values):

	#plot vars set up
	genebody = config['region'] == 'genebody'
	plottype = config['plottype']

	if genebody:
		size = 3*flank
		s = flank/10000
	else:
		size = 2*flank+1
		s = flank/10000
	
	nticks = config['nticks']
	if genebody:
		ticks = ['-'+str(i*flank//nticks) for i in reversed(range(1,nticks+1))]+['TSS']+['%.2f' % float(i/(nticks)) for i in range(1,nticks)]+['TSE']+['+' + str(i*flank//nticks) for i in range(1,nticks+1)]
		if nticks == 0:
			tickvals = [flank,2*flank]
		else:
			tickvals = [i*flank//nticks for i in range(0,nticks)]+[flank+i*flank//(nticks) for i in range(0,nticks)]+[2*flank+i*flank//nticks for i in range(0,nticks+1)]
	else:
		ticks = ['-'+str(i*flank//nticks) for i in reversed(range(1,nticks+1))]+[region.upper()]+['+' + str(i*flank//nticks) for i in range(1,nticks+1)]
		if nticks == 0:
			tickvals = [flank]
		else:
			tickvals = np.linspace(0,size,2*nticks+1)

	if config['smooth']:
		if config['smooth'].lower() == 'true':
			smooth = flank/10000
		elif config['smooth'].lower() == 'false':
			config['smooth'] = False
		else:
			smooth = float(config['smooth'])


	#plot avgprof
	if plottype in ['avgprof','both']:

		print('calculating mean...')

		if config['smooth']:
			spl = interpolate.UnivariateSpline(np.linspace(0,size,size),np.mean(values, axis=0,dtype = 'float16'),s = smooth)
			avgprof = spl(np.linspace(0,size,size))

		else:
			avgprof = np.mean(values, axis=0,dtype = 'float16')

		print('plotting...')
		with sb.axes_style("darkgrid"):
			plt.plot(np.linspace(0,size,size),avgprof)

			plt.xticks(tickvals,ticks)
			plt.xlim((0,size))
			plt.ylim((avgprof.min()-0.05,avgprof.max()*1.05))
			if config['avgtitle']:
				plt.title(config['avgtitle'])

			if config['avgfile']:
				plt.savefig(config['avgfile'])
			else:
				plt.show()
	#plot heatmap
	if plottype in ['heatmap','both']:

		if config['sort']:

			print('sorting...')
			b = np.sum(values,axis = 1) * -1
			indx = b.argsort()
			values = np.take(values,indx,axis=0)

		if plottype == 'both':
			plt.figure()

		print('plotting...')
		with sb.axes_style("ticks"):
			if config['hnorm'] == 'lin':
				plt.imshow(values,aspect = 'auto',cmap = config['cmap'],norm = colors.Normalize(vmin=values.min(),vmax=values.max()*.8))
			elif config['hnorm'] == 'log':
				plt.imshow(values,aspect = 'auto',cmap = config['cmap'],norm =colors.LogNorm())

			plt.xticks(tickvals,ticks)
			plt.xlim((0,size))

			if config['cbar']:
				plt.colorbar()
			if config['hmtitle']:
				plt.title(config['hmtitle'])

			if config['hmfile']:
				plt.savefig(config['hmfile'])
			else:
				plt.show()
			
def readScore():
	f = open(config['gfile'],'r')
	result = []

	if config['gomit']:
		gOmit = config['gomit']
	else:
		gOmit = ()

	if config['chrmomit']:
		cOmit = config['chrmomit']
	else:
		cOmit = ()

	if config['chrmonly']:
		cOnly = config['chrmonly']
	else:
		cOnly = None

	for line in f:
		l = line.split('_')
		gName = l[1]
		append = True
		for x in gOmit:
			if gName.startswith(x):
				append = False
				break
		l = l[-1].split(':')
		chrm = l[0][:-1]
		if cOnly:
			for x in cOnly:
				if chrm not in cOnly:
					append = False
		if append and chrm not in cOmit:
			l = l[1].split('\t')
			score = int(l[1])
			l = l[0].split('-')
			result.append((chrm,int(l[0]),int(l[1]),gName,score))
	f.close()

	if config['nbest']:
		n = config['nbest']
		if n < len(result):
			result = sorted(result,key = lambda t: t[-1],reverse = True)
			result = result[:n]
	elif config['scorerange']:
		r = config['scorerange']
		r = (int(r[0]),int(r[1]))
		result = sorted(result,key = lambda t: t[-1])
		scores = [t[-1] for t in result]
		i = bisect_left(scores,r[0])
		j = bisect_right(scores,r[1],lo = i)
		result = result[i:j]

	return result

def readBed():
	
	if config['gindx']:
		l = config['gindx']
		indx = {'chrm':int(l[0]),'start':int(l[1]),'end':int(l[2]),'name':int(l[3]),'strand':int([4])}
	else:
		indx = {'chrm':0,'start':1,'end':2,'name':3,'strand':-1}

	f = open(config['gfile'],'r')
	result = []

	if config['gomit']:
		gOmit = config['gomit']
	else:
		gOmit = ()
	if config['chrmomit']:
		cOmit = config['chrmomit']
	else:
		cOmit = ()

	ofirst = config['ofirst']

	for line in f:
		line = line.split()
		if line[indx['chrm']] not in cOmit:
			append = True
			for g in gOmit:
				if line[indx['name']].startswith(g):
					append = False
					break

			if ofirst and len(result)>0:
				if result[-1][1] == line[indx['name']]:
					append = False
			if append:
				if line[-1] == '+':
					result.append((line[indx['chrm']],int(line[indx['start']]),int(line[indx['end']]),line[indx['name']]))
				else:
					result.append((line[indx['chrm']],int(line[indx['end']]),int(line[indx['start']]),line[indx['name']]))
	f.close()
	return result

readConfig()
region = config['region']
flank = config['flank']
TSE = region =='tse'

if config['replot']:
	print('loading file')
	plot(np.load(config['replot']))
elif config['gfile'] and config['infile']:
	s = loadGenes()
	t = query(config['infile'],s)
	if config['matfile']:
		np.save(config['matfile'],t)
	print('array generated')
	plot(t)
	print('Done')
else:
	print("Input files not provided")
