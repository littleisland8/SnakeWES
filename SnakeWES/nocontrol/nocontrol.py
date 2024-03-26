#!/usr/bin/python3 env

import subprocess 
import shlex 
import yaml

class c():

	'''
	Container. This stores argparser parameters. Used to pass multiple parameters at once.
	'''

	OneSample=False
	threads=''

def loadConfig():

	with open("SnakeWES/config/conf.yaml","r") as fin:

		data=yaml.load(fin,Loader=yaml.FullLoader)

	return data


def run(parser,args):

	c.OneSample = args.OneSample
	c.threads= args.threads

	configdict=loadConfig()
	configdict["run_mode"]="nocontrol"

	with open("SnakeWES/config/conf.yaml","w") as fout:

		yaml.dump(configdict, fout, default_flow_style=False)

	if c.OneSample:

		command="snakemake NocontrolsOneSample --cores " + str(c.threads) + " --use-conda"

	else:

		command="snakemake Nocontrols --cores " + str(c.threads) + " --use-conda"

	subprocess.call(shlex.split(command))