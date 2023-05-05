#!/usr/bin/env python3
# coding: utf-8
"""
.. module:: run_BUSCO
   :synopsis:
.. versionadded:: 3.0.0
.. versionchanged:: 5.4.0

BUSCO - Benchmarking Universal Single-Copy Orthologs.
This is the BUSCO main script.

To get help, ``busco -h``. See also the user guide.

And visit our website `<http://busco.ezlab.org/>`_

Copyright (c) 2016-2023, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""


import argparse
from argparse import RawTextHelpFormatter
from libs.BUSCO.src.busco._version import __version__
from libs.BUSCO.src.busco.BuscoRunner import AnalysisRunner, BatchRunner, SingleRunner
from libs.BUSCO.src.busco.Exceptions import BatchFatalError, BuscoError
from libs.BUSCO.src.busco.BuscoLogger import BuscoLogger
from libs.BUSCO.src.busco.BuscoLogger import LogDecorator as log
from libs.BUSCO.src.busco.ConfigManager import BuscoConfigManager
from libs.BUSCO.src.busco.Actions import (
    ListLineagesAction,
    CleanHelpAction,
    CleanVersionAction,
    DirectDownload,
)
from libs.BUSCO.src.busco.ConfigManager import BuscoConfigMain

import sys
import time


logger = BuscoLogger.get_logger(__name__)


@log(
    "***** Start a BUSCO v{} analysis, current time: {} *****".format(
        __version__, time.strftime("%m/%d/%Y %H:%M:%S")
    ),
    logger,
)
class BuscoMaster:
    def __init__(self, params):
        self.params = params
        self.config_manager = BuscoConfigManager(self.params)
        self.config = self.config_manager.config_main

    def load_config(self):
        """
        Load a busco config file that will figure out all the params from all sources
        i.e. provided config file, dataset cfg, and user args
        """
        self.config_manager.load_busco_config_main()
        self.config = self.config_manager.config_main

    def check_batch_mode(self):
        return self.config.getboolean("busco_run", "batch_mode")

    def run(self):
        try:
            self.load_config()
            runner = (
                BatchRunner(self.config_manager)
                if self.check_batch_mode()
                else SingleRunner(self.config_manager)
            )
            out = runner.run()

        except BuscoError as be:
            SingleRunner.log_error(be)
            raise SystemExit(1)

        except BatchFatalError as bfe:
            SingleRunner.log_error(bfe)
            raise SystemExit(1)

        finally:
            try:
                AnalysisRunner.move_log_file(self.config)
            except:
                pass
        return out


@log("Command line: {}".format(" ".join(sys.argv[:])), logger, debug=True)
def parameter(assembly, outputdir, lineage, auto_l):
    param = {'in': assembly, 'out': outputdir, 'mode': 'genome',
     'lineage_dataset': lineage, 'use_augustus': False, 'augustus_parameters': None, 'augustus_species': None,
     'auto-lineage': auto_l, 'auto-lineage-euk': False, 'auto-lineage-prok': False, 'cpu': 8, 'config_file': None,
     'contig_break': None, 'datasets_version': None, 'download': '==SUPPRESS==', 'download_base_url': None,
     'download_path': None, 'evalue': None, 'force': False, 'help': '==SUPPRESS==', 'limit': None,
     'list_datasets': '==SUPPRESS==', 'long': False, 'metaeuk_parameters': None, 'metaeuk_rerun_parameters': None,
     'offline': False, 'out_path': None, 'quiet': False, 'restart': False, 'scaffold_composition': False, 'tar': False,
     'update-data': False, 'version': '==SUPPRESS=='}
    return param


def run_busco(assembly, outputdir, lineage, auto_l):
    """
    This function runs a BUSCO analysis according to the provided parameters.
    See the help for more details:
    ``busco -h``
    :raises SystemExit: if any errors occur
    """
    params = parameter(assembly, outputdir, lineage, auto_l)
    busco_run = BuscoMaster(params)
    return busco_run.run()


