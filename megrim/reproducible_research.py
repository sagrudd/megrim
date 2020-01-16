
import sys
import logging
import re
import os
import pandas as pd
import platform
import locale
import math


class PythonModule:

    internal_modules = ["marshal", "sys"]

    def __init__(self, name, version):
        self.name = name
        head, sep, tail = str(version).partition(' ')
        self.version = head

    def get_name(self):
        return self.name

    def get_version(self):
        return self.version

    def is_internal_module(self):
        return self.name in self.internal_modules


class SessionInfo:

    def __init__(self, include_internal=False):
        self.include_internal = include_internal
        self.modules_found = []
        self.parse_modules()

    def parse_modules(self, include_na=False, return_list=False):

        for module in sys.modules:
            logging.debug(module)
            if '.' not in module and not module.startswith('_'):
                if hasattr(sys.modules[module], '__version__'):
                    self.module_append(PythonModule(module, sys.modules[module].__version__))
                elif hasattr(sys.modules[module], 'VERSION'):
                    self.module_append(PythonModule(module, sys.modules[module].VERSION))
                elif hasattr(sys.modules[module], 'version'):
                    if callable(sys.modules[module].version):
                        self.module_append(PythonModule(module, sys.modules[module].version()))
                    else:
                        self.module_append(PythonModule(module, sys.modules[module].version))

    def module_append(self, module):
        if not module.is_internal_module() or self.include_internal:
            self.modules_found.append(module)

    def list_modules(self):
        module_names = []
        for module in self.modules_found:
            module_names.append(module.get_name())
        module_names = sorted(module_names)
        return module_names

    def get_module(self, module_name):
        for module in self.modules_found:
            if module.get_name() == module_name:
                return module
        return None

    def get_version_list(self, module_names):
        module_versions = []
        for module_name in module_names:
            module_versions.append(self.get_module(module_name).get_version())
        return module_versions

    def get_module_df(self):
        module_names = self.list_modules()
        d = {'module_name': module_names, 'module_version': self.get_version_list(module_names)}
        return pd.DataFrame(data=d)

    def list_vmodules(self, pad=0):
        module_names = self.list_modules()
        module_vnames = []
        for module_name in module_names:
            module_vnames.append("%s {%s}%s" % (module_name, self.get_module(module_name).get_version(), ' '*pad))
        return module_vnames

    def get_simple_df(self, cols=3, pad=0):
        module_vnames = self.list_vmodules(pad=pad)
        items = len(module_vnames)
        delta = items - math.floor(items / cols) * cols
        if delta > 0:
            add_items = cols - delta
            for i in range(add_items):
                module_vnames.append('')
        df = pd.DataFrame(data=pd.Series(module_vnames).values.reshape((-1, cols), order="F"))
        return df

    def get_simple_df_str(self, cols=4, prefix="   ", pad=3):
        df = self.get_simple_df(cols=cols, pad=pad)
        rows = df.to_string().split("\n")[1:]
        stripc = 0
        if 0 <= len(df.index) < 10:
            stripc = 1
        if 0 <= len(df.index) < 10:
            stripc = 2
        newrows = []
        for row in rows:
            newrows.append(prefix+row[stripc:])
        return "attached modules:\n" + "\n".join(newrows)

    def get_python_version(self):
        return sys.version.partition('\n')[0].strip()

    def get_python_version_str(self):
        return "Python version: %s" % self.get_python_version()

    def get_platform(self):
        return "-".join([platform.machine(), platform.system()])


    def get_platform_str(self):
        return "Platform: %s" % self.get_platform()

    def get_os(self):
        return platform.platform()

    def get_os_str(self):
        return "Running under: %s" % self.get_os()

    def session_info(self):
        return "\n".join([self.get_python_version_str(), self.get_platform_str(), self.get_os_str(), '',
                          self.get_locale_str(), self.get_simple_df_str()])

    def get_locale_str(self):
        return "locale: %s" % str(locale.getlocale())

    def __str__(self):
        return self.session_info()


print(SessionInfo())

