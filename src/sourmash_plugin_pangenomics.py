"""pangenomics plugin"""

epilog="""
See https://github.com/dib-lab/sourmash_plugin_pangenomics for more examples.

Need help? Have questions? Ask at http://github.com/sourmash-bio/sourmash/issues!
"""

import argparse
import sourmash

from sourmash.logging import debug_literal
from sourmash.plugins import CommandLinePlugin

###

#
# CLI plugins - supports 'sourmash scripts <commands>'
#

class Command_CreateDB(CommandLinePlugin):
    command = 'pangenomics_createdb'             # 'scripts <command>'
    description = "@CTB"       # output with -h
    usage = "@CTB"               # output with no args/bad args as well as -h
    epilog = epilog             # output with -h
    formatter_class = argparse.RawTextHelpFormatter # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        # add argparse arguments here.
        debug_literal('RUNNING cmd_xyz.__init__')

    def main(self, args):
        # code that we actually run.
        super().main(args)
        print('RUNNING cmd', self, args)


class Command_RankTable(CommandLinePlugin):
    command = 'pangenomics_ranktable'             # 'scripts <command>'
    description = "@CTB"       # output with -h
    usage = "@CTB"               # output with no args/bad args as well as -h
    epilog = epilog             # output with -h
    formatter_class = argparse.RawTextHelpFormatter # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        # add argparse arguments here.
        debug_literal('RUNNING cmd_xyz.__init__')

    def main(self, args):
        # code that we actually run.
        super().main(args)
        print('RUNNING cmd', self, args)


class Command_Classify(CommandLinePlugin):
    command = 'pangenomics_classify'             # 'scripts <command>'
    description = "@CTB"       # output with -h
    usage = "@CTB"               # output with no args/bad args as well as -h
    epilog = epilog             # output with -h
    formatter_class = argparse.RawTextHelpFormatter # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        # add argparse arguments here.
        debug_literal('RUNNING cmd_xyz.__init__')

    def main(self, args):
        # code that we actually run.
        super().main(args)
        print('RUNNING cmd', self, args)
