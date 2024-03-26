#!/usr/bin/python3 env

import sys
import argparse
from argparse import HelpFormatter

#from Test import __version__

def main():

	parser = argparse.ArgumentParser(prog='SnakeWES', description='''Snakemake pipeline for WES analysis''', epilog='''This program was developed by Simone Romagnoli and Alessio Enderti''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='PoN, Paired, Nocontrol')

	## PoN ##

	parser_PoN = subparsers.add_parser('PoN', help='Perform WES analysis considering a panel of normal samples')

	required = parser_PoN.add_argument_group('Required I/O arguments')
	required.add_argument('--threads', help='number of cores to use for Snakemake pipeline', required=True, metavar='', type=int, default=1)
	
	additionals=parser_PoN.add_argument_group('Additional arguments')
	additionals.add_argument('--OneSample', help='Perform PoN analysis considering 1 tumor sample', action='store_true', default=False)	

	parser_PoN.set_defaults(func=run_subtool)

	## Paired ##
	
	parser_paired = subparsers.add_parser('Paired', help='Perform WES analysis considering paired germline-tumor samples')

	required = parser_paired.add_argument_group('Required I/O arguments')
	required.add_argument('--threads', help='number of cores to use for Snakemake pipeline', required=True, metavar='', type=int, default=1)
	
	additionals=parser_paired.add_argument_group('Additional arguments')
	additionals.add_argument('--OneSample', help='Perform Paired analysis considering 1 tumor sample', action='store_true', default=False)	

	parser_paired.set_defaults(func=run_subtool)

	## Nocontrol ##
	
	parser_nocontrol = subparsers.add_parser('Nocontrol', help='Perform WES analysis considering only tumor samples')

	required = parser_nocontrol.add_argument_group('Required I/O arguments')
	required.add_argument('--threads', help='number of cores to use for Snakemake pipeline', required=True, metavar='', type=int, default=1)
	
	additionals=parser_nocontrol.add_argument_group('Additional arguments')
	additionals.add_argument('--OneSample', help='Perform No control analysis considering 1 tumor sample', action='store_true', default=False)

	parser_nocontrol.set_defaults(func=run_subtool)

	if len(sys.argv)==1:
    	
		parser.print_help(sys.stderr)
		sys.exit(1)

	#case-insensitive submodules
	
	if sys.argv[1].lower() == 'pon':

		sys.argv[1] = 'PoN'

	elif sys.argv[1].lower() == 'paired':

		sys.argv[1] = 'Paired'

	elif sys.argv[1].lower() == 'nocontrol':

		sys.argv[1] = 'Nocontrol'

	args = parser.parse_args()
	args.func(parser, args)


class CustomFormat(HelpFormatter):

	'''
	Customize how help is diplayed
	'''

	def _format_action_invocation(self, action):

		if not action.option_strings:

			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			
			return metavar

		else:

			parts = []

			if action.nargs == 0:

				parts.extend(action.option_strings)

			else:

				default = self._get_default_metavar_for_optional(action)
				args_string = self._format_args(action, default)
				
				for option_string in action.option_strings:

					parts.append(option_string)

				return '%s %s' % (', '.join(parts), args_string)

			return ', '.join(parts)

	def _get_default_metavar_for_optional(self, action):

		return action.dest.upper()


def run_subtool(parser, args):


	if args.command == 'PoN': 

		from .PoN import PoN as submodule
	
	elif args.command == 'Paired': 

		from .Paired import Paired as submodule

	elif args.command == 'Nocontrol':

		from .Nocontrol import Nocontrol as submodule

	else:

		parser.print_help()

	submodule.run(parser,args)


if __name__ =='__main__':

	main()
