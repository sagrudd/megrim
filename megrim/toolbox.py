"""
Develop functionality for megrim commandline plugins.

This module aims to provide some functionality within the broader megrim
package to enable command line access to functionality through an extensible
plugin framework.

Created on Thu Feb 20 13:22:55 2020

@author: srudd
"""

import inspect
import pkgutil
from megrim.environment import MegrimPlugin
import argparse
from importlib import reload
import logging
import tempfile
import multiprocessing
import warnings


class MegrimToolBox:
    """The toolbox for providing generic functionality."""

    def __init__(self, plugin_package, target="plugins"):
        self.plugin_package = plugin_package
        self.target = target
        self.plugins = []
        self.reload_plugins()

    def reload_plugins(self):
        """
        Reload the available plugins within the working environment.

        Returns
        -------
        None.

        """
        logging.debug(
            f"looking for megrim plugins in package {self.plugin_package}")
        if self.target is None:
            self.walk_package(self.plugin_package)
        else:
            self.walk_package(f"{self.plugin_package}.{self.target}")

    def walk_package(self, package):
        """
        Attempt to load an identified plugin module.

        Parameters
        ----------
        package: str
            The package that should be loaded.

        Returns
        -------
        None.

        """
        # print(f'Walk-package ... {package}')
        imported_package = __import__(package, fromlist=[""])
        for _, pluginname, ispkg in pkgutil.iter_modules(
                imported_package.__path__, imported_package.__name__ + "."):
            if not ispkg:
                # print(f"checking {pluginname}")
                plugin_module = __import__(pluginname, fromlist=[""])
                # print(plugin_module)
                clsmembers = inspect.getmembers(plugin_module, inspect.isclass)
                for (_, c) in clsmembers:
                    # print(f"{_}\t\t{c}")
                    if issubclass(c, MegrimPlugin) and (c is not MegrimPlugin):
                        logging.debug(
                            f'\timporting plugin class: {c.__name__}')
                        self.plugins.append(c())

    def list(self):
        """
        List the available plugins.

        Returns
        -------
        list
            The available plugins in a comma-delimitted string format.

        """
        code_words = []
        for plugin in self.plugins:
            code_words.append(plugin.tool)
        return ", ".join(code_words)

    def execute(self, args):
        """
        Run the execute method of the specified plugin.

        Parameters
        ----------
        args: TYPE
            DESCRIPTION.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        try:
            fubar = True
            for plugin in self.plugins:
                if plugin.tool == args.method:
                    fubar = False
                    plugin.execute(args)
            if fubar:
                raise ValueError(
                    f"the requested method [{args.method}] is not known")
        except ValueError as e:
            print("Exception!", e)

    def arg_params(self, subparsers, parent_parser):
        """
        Pull configuration parameters defined within the available plugins.

        Parameters
        ----------
        subparsers : TYPE
            DESCRIPTION.
        parent_parser : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        logging.debug("merging annotations ...")
        for plugin in self.plugins:
            plugin.arg_params(subparsers, parent_parser)


def main():
    """
    Provide a main method to orchestrate plugin framework and run stuff.

    Returns
    -------
    None.

    """
    warnings.simplefilter(action='ignore', category=FutureWarning)
    reload(logging)
    logging.basicConfig(
        format='%(asctime)s %(levelname)s:%(message)s',
        level=logging.INFO, datefmt='%I:%M:%S')

    megrim_plugins = MegrimToolBox("megrim")
    parser = argparse.ArgumentParser()
    # parser.add_argument("method", help=f"Define the megrim method to run]")

    subparsers = parser.add_subparsers(title='subcommand help')
    subparsers.required = True
    subparsers.dest = 'method'

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        '--cache', metavar="/tmp", action='store',
        help='Path to location for storing cached and temporary files.',
        dest="cache", default=tempfile.gettempdir())
    parent_parser.add_argument(
        '--debug', action='store_true', dest="debug",
        help='Display debug information and messages', default=False)
    parent_parser.add_argument(
        '--threads', action='store', dest="threads",
        help=f'The number of threads to use for parallel steps of a workflow '
        '(default={max(1, int(multiprocessing.cpu_count()/2))})',
        default=max(1, int(multiprocessing.cpu_count()/2)), type=int)

    megrim_plugins.arg_params(subparsers, parent_parser)
    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    megrim_plugins.execute(args)
