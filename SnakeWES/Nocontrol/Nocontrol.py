#!/usr/bin/python3 env

import subprocess 
import shlex 
import yaml

class c():

	'''
	Container. This stores argparser parameters. Used to pass multiple parameters at once.
	'''

	OneSample=False
	dr=False
	threads=''

def loadConfig():

	with open("SnakeWES/config/conf.yaml","r") as fin:

		data=yaml.load(fin,Loader=yaml.FullLoader)

	return data


def run(parser,args):

	c.OneSample = args.OneSample
	c.threads= args.threads
	c.dr=args.dr

	configdict=loadConfig()
	configdict["run_mode"]="Nocontrol"

	with open("SnakeWES/config/conf.yaml","w") as fout:

		yaml.dump(configdict, fout, default_flow_style=False)

	if c.OneSample:

		if c.dr:

			command="snakemake NocontrolsOneSample --cores " + str(c.threads) + " --use-conda -np"

		else:

			command="snakemake NocontrolsOneSample --cores " + str(c.threads) + " --use-conda"	

	else:

		if c.dr:

			command="snakemake Nocontrols --cores " + str(c.threads) + " --use-conda -np"

		else:

			command="snakemake Nocontrols --cores " + str(c.threads) + " --use-conda"

	subprocess.call(shlex.split(command))