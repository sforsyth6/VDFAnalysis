import os
from argparse import ArgumentParser

import time

from astropy.table import Table
from astropy.io.ascii import (write,read)

from gwpy.time import to_gps
from gwpy.segments import (DataQualityFlag, DataQualityDict)
from gwpy.segments.segments import SegmentList

from gwpy.table.lsctables import SnglBurstTable
from gwpy.table import EventTable


from gwsumm.utils import mkdir

from gwvet.segments import get_known_flags
from gwvet.metric.metrics import (deadtime, efficiency)
from gwvet.metric.__init__ import evaluate
from gwvet.metric.registry import get_metric
from urlparse import urlparse
from urllib2 import urlopen

start_time = time.time()

parser = ArgumentParser(description=__doc__)
parser.add_argument('-v', '--verbose', action='store_true',
                    help='print verbose output')

# required argument
parser.add_argument('veto-definer-file', help='path to veto definer file')
parser.add_argument('gps-start-time', type=to_gps,
                    help='GPS start time/date of analysis')
parser.add_argument('gps-end-time', type=to_gps,
                    help='GPS end time/date of analysis')

analargs = parser.add_argument_group('Analysis options')
analargs.add_argument('-i', '--ifo',
                      help='prefix of IFO to study, default: %(default)s')
analargs.add_argument('-a', '--analysis', default = 'omciron',
		      help='either omicron or cWb')
analargs.add_argument('-t', '--trigFile', default = '/home/steven.forsyth/public_html/analysis/flagAnalysis/cWb/H1/O1/O1triggers/constrained/EVENTS.txt',
		      help='location of your triggers text file')
analargs.add_argument('-o', '--output-directory', default=os.curdir,
                      type=os.path.abspath,
                      help='output directory path, default: %(default)s, '
                           'this path should be web-viewable')

url = 'https://segments.ligo.org'
args = parser.parse_args()
start = getattr(args, 'gps-start-time')
end = getattr(args, 'gps-end-time')
ifo = args.ifo
vetofile = getattr(args, 'veto-definer-file')
metrics = ['efficiency','efficiency/deadtime','Efficiency/Deadtime | SNR>=8', 'Efficiency/Deadtime | SNR>=20', 'Efficiency/Deadtime | SNR>=100']
trigPath = args.trigFile
analysis = args.analysis

#create directories to hold stuff
tag = '%d-%d-%s' % (start.seconds, end.seconds,analysis)
outdir = os.path.abspath(os.path.join(args.output_directory, tag))
mkdir(outdir)
os.chdir(outdir)

ALLSEGMENTS = DataQualityDict()

# -- get analysis flags ----------------------
allFlags = get_known_flags(start,end,url, ifo=ifo, badonly=None)
for flag in allFlags[:]:
    if 'ANALYSIS_READY' not in flag:
	new = DataQualityFlag.query(flag,start,end) #try get_segments gwvet.segments
    	ALLSEGMENTS[new.name] = new
    else:
       allFlags.remove(flag)

#Download VDFs
if urlparse(vetofile).netloc:
    tmp = urlopen(vetofile)
    vetofile = os.path.abspath(os.path.basename(vetofile))
    with open(vetofile, 'w') as f:
        f.write(tmp.read())
    print('Downloaded veto definer file')
vdf = DataQualityDict.from_veto_definer_file(
    vetofile, format='ligolw', start=start, end=end, ifo=ifo)
vdf.populate()
print('Read %d flags from veto definer' % len(vdf.keys()))

categories = [1,2,3,4]
vdfCats = dict((c, DataQualityDict()) for c in categories)
for name, flag in vdf.iteritems():
    try:
        vdfCats[flag.category][name] = flag
    except KeyError:
        pass

#generate the triggers needed to run vet
if analysis == 'omicron':
	triggers = SnglBurstTable.fetch('%s:GDS-CALIB_STRAIN' %ifo, 'omicron', start, end)
else:
	triggers = SnglBurstTable.read(trigPath,ifo=ifo,format='cwb-ascii')

#convert the metrics to a standard notation
newMetric = []
for metric in metrics:
	newMetric.append(get_metric(metric))

order = ['1','1+4','1+2+4', '1+2+3+4']
for categories in order:
	category = categories.split('+')
	activeVdf = DataQualityFlag('VDFs')
	for cat in category:
		for flag in vdfCats[int(cat)]:
			activeVdf.known += vdf[flag].known
			activeVdf.active += vdf[flag].active

	vetoTrigs = triggers.veto(activeVdf.active)
	withoutVdf = {'Flags': [], 'Deadtime':[], 'Efficiency':[],'Eff/ded': [], 'E/D SNR >= 8':[], 'E/D SNR >= 20':[], 'E/D SNR >= 100':[]}

	for flags in ALLSEGMENTS.keys():
		segment = ALLSEGMENTS[flags]
		values = evaluate(segment,triggers,newMetric)
		
		valuesWo = evaluate(segment,vetoTrigs,newMetric)
		withoutVdf['Flags'].append(flags)
		withoutVdf['Deadtime'].append(deadtime(segment))
		withoutVdf['Efficiency'].append(valuesWo[0]*100)
		withoutVdf['Eff/ded'].append(valuesWo[1])
		withoutVdf['E/D SNR >= 8'].append(valuesWo[2])
		withoutVdf['E/D SNR >= 20'].append(valuesWo[3])
		withoutVdf['E/D SNR >= 100'].append(valuesWo[4])	

	#create a table of the values
	withoutVdfTable = Table(withoutVdf,names = ('Flags','Deadtime','Efficiency', 'Eff/ded', 'E/D SNR >= 8', 'E/D SNR >= 20', 'E/D SNR >= 100') )

	write(withoutVdfTable,'cat%s.txt' %(categories), format='fixed_width_two_line')
