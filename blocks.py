import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('blocks')
parser.add_argument('ifo')

args = parser.parse_args()
ifo = getattr(args,'ifo')
blocks = int(getattr(args,'blocks'))

start = 1126051217
start_time = start
end_segment = [1127271617,1128299417,1129383017,1130754617,1132104617,1133173817,1134450017, 1135652417, 1137254417]

i=0
cmd = []
while i <= blocks:
	if i == blocks:
		start_time = start
		end_time = end_segment[i-1]
	else:
		end_time = end_segment[i]
	
	cmd.append('python analysis.py -i %s -a cWb https://code.pycbc.phy.syr.edu/detchar/veto-definitions/download/master/burst/O1/H1L1-HOFT_C02_O1_BURST.xml %i %i' %(ifo,start_time,end_time))
	start_time = end_time
	i +=1

i=1
code = cmd[0]
while i <= blocks:
	code += ' & %s' %cmd[i]
	i += 1
os.system(code)
