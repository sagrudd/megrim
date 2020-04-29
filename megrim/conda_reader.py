import logging
from megrim.environment import Flounder
import argparse
import yaml
import os
import re
import urllib.request
import ssl
from urllib.parse import urlparse
import sys


class CondaGit(Flounder):

    def __init__(self, args):
        Flounder.__init__(self)

        if args is not None:
            if isinstance(args, argparse.Namespace):
                self.argparse(args)
            if isinstance(args, dict):
                self.dictparse(args)

        self.conda_yaml = None

    def lookup(self):
        self.extract_conda_manifest()

    def extract_conda_manifest(self):
        logging.info(f"looking up conda entity [{self.args.target}] from [{self.args.bioconda}] ...")

        # read the YAML file ...
        src = os.path.join(self.args.bioconda, "recipes", self.args.target, "meta.yaml")
        logging.info(f"looking for {src}")
        # read yaml lines as text ... strip out the {{ set whatevers
        with open(src) as file:
            lines = file.readlines()
        # and parse the lines ...
        substitutions = {}
        text = ""
        for line in lines:
            if re.search("^{% set", line):
                field, value = extract_set_value(line)
                substitutions[field] = value
            else:
                text = text + correct_set_line(line, substitutions)
        self.conda_yaml = yaml.load(text, Loader=yaml.BaseLoader)
        logging.info(self.conda_yaml)

    def get_build_lines(self):
        build = os.path.join(self.args.bioconda, "recipes", self.args.target, "build.sh")
        logging.info(f"looking for {build}")
        with open(build) as file:
            build_lines = file.readlines()

        processed_lines = []
        # we should collapse the multi-line \ entries

        entry = ""
        for line in build_lines:
            line = line[:-1]
            if line.endswith("\\"):
                entry += line[:-1].strip() + " "
            elif len(entry) > 0:
                processed_lines.append(entry)
                entry = ""
            else:
                processed_lines.append(line)

        return processed_lines


def extract_set_value(line):
    line = re.sub("({% set|%}|\")", "", line.strip()).strip().split("=")
    logging.info(f"extracted set value == {line[0].strip()} ... {line[1].strip()}")
    return line[0].strip(), line[1].strip()


def correct_set_line(line, substitutions):
    for key in substitutions.keys():
        line = re.sub("".join(["{{ ", key, " }}"]), substitutions[key], line)
    line = re.sub("{{" , "\"{{", line)
    line = re.sub("}}", "}}\"", line)
    return line



class RpmHandler(Flounder):

    def __init__(self, args, conda):
        Flounder.__init__(self)

        if args is not None:
            if isinstance(args, argparse.Namespace):
                self.argparse(args)
            if isinstance(args, dict):
                self.dictparse(args)
        self.conda = conda
        self.has_bin = False
        self.defined_vals = {"{CPU_COUNT}": "4", "SRC_DIR": "%{buildroot}"}

    def prepare_manifest(self):

        if not os.path.exists(self.args.rpm):
            os.makedirs(self.args.rpm)
            os.makedirs(os.path.join(self.args.rpm, "BUILD"))
            os.makedirs(os.path.join(self.args.rpm, "RPMS"))
            os.makedirs(os.path.join(self.args.rpm, "SOURCES"))
            os.makedirs(os.path.join(self.args.rpm, "SPECS"))
            os.makedirs(os.path.join(self.args.rpm, "SRPMS"))

        package_name = self.conda.conda_yaml["package"]["name"]
        package_version = self.conda.conda_yaml["package"]["version"]
        package_release = self.conda.conda_yaml["build"]["number"]
        spec_id = f"{package_name}-{package_version}.spec"
        spec = os.path.join(self.args.rpm, "SPECS", spec_id)

        if os.path.exists(spec):
            logging.info(f"RPM manifest [{spec_id}] already exists")
            logging.info(f"             release version at [{package_release}]")
            if not self.args.force:
                logging.info("Skipping this file ...")
                return
            else:
                logging.info("incrementing release version for RPM updates")
                spec_pars = SpecParser(spec)
                existing_release = spec_pars.get_spec_release()
                if int(existing_release) >= int(package_release):
                    package_release = int(existing_release) + 1
                logging.info(f"RPM SPEC being updated with release version [{package_release}]")

        logging.info(f"preparing RPM manifest [{spec}]")
        with open(spec, "w") as file:
            print("%define debug_package %{nil}", file=file)
            print("\nName: "+package_name, file=file)
            print("Version: "+package_version, file=file)
            print("Release: "+str(package_release)+"%{?dist}", file=file)
            print("Summary: "+self.conda.conda_yaml["about"]["summary"], file=file)
            print("Group: Applications / Bioinformatics", file=file)
            print("License: "+self.conda.conda_yaml["about"]["license"], file=file)
            print("URL: "+self.conda.conda_yaml["about"]["home"], file=file)
            print("Source0: "+self.conda.conda_yaml["source"]["url"], file=file)

            build_requires = []
            for build in self.conda.conda_yaml["requirements"]["build"]:
                if build in conda2rpm_mapping.keys():
                    build = conda2rpm_mapping[build]
                if build not in build_requires:
                    build_requires.append(build)

            for build in self.conda.conda_yaml["requirements"]["host"]:
                if build in conda2rpm_mapping.keys():
                    build = conda2rpm_mapping[build]
                if build not in build_requires:
                    build_requires.append(build)

            run_requires = ["environment-modules"]
            for build in self.conda.conda_yaml["requirements"]["run"]:
                if build in conda2rpm_mapping.keys():
                    build = conda2rpm_mapping[build]
                if build not in run_requires:
                    run_requires.append(build)

            print("\nBuildRequires: "+" ".join(build_requires), file=file)
            print("Requires: "+" ".join(run_requires), file=file)

            print("\n%description\n"+self.conda.conda_yaml["about"]["summary"], file=file)

            print("%define _bindir /bin", file=file)
            print("%define _libdir /lib", file=file)
            print("%define _mandir /man", file=file)
            print("%define _datarootdir /", file=file)

            print("\n%prep", file=file)

            print("\n%setup", file=file)

            # we need some configure handling to be added here?
            # print("\n%configure", file=file)
            self.extract_configuration_cmds(file)

            print("\n%build", file=file)
            self.extract_build_cmds(file)

            print("\n%install", file=file)
            self.extract_install_cmds(file)

            print("\n%files", file=file)
            self.dump_file_content(file)

            print("\n%clean", file=file)
            self.cleanup(file)

            print("\n%changelog", file=file)


        self.download_source()

    def extract_configuration_cmds(self, fh):
        #print("=====")
        build_lines = self.conda.get_build_lines()
        for line in build_lines:
            line = line.strip()
            line = self.parseprefix(line)
            if re.search("^[^\\s]+=", line):
                key, value = line.split("=", maxsplit=1)
                if key in self.defined_vals.keys() and re.search("\\$"+key, value):
                    value = re.sub("(\\$"+key+"|\"\")", "", self.defined_vals[key] + value)

                self.defined_vals[key] = value
                logging.info(f"internal value [{key}] --> [{value}]")


    def extract_build_cmds(self, fh):
        build_lines = self.conda.get_build_lines()
        for line in build_lines:
            line = line.strip()
            line = self.parseprefix(line)

            for key in self.defined_vals.keys():
                if re.search("\\$"+key, line):
                    line = re.sub("\\$"+key, self.defined_vals[key], line)

            # print(line)
            if re.search("^./configure", line):
                print(line, file=fh)
            if re.search("^make", line) and not re.search("^make install", line):
                print(line, file=fh)
            if re.search("^curl", line):
                print(line, file=fh)




    def extract_install_cmds(self, fh):
        build_lines = self.conda.get_build_lines()
        for line in build_lines:
            line = line.strip()
            line = self.parseprefix(line)

            for key in self.defined_vals.keys():
                if re.search("\\$"+key, line):
                    line = re.sub("\\$"+key, self.defined_vals[key], line)

            print(line)
            if re.search("^mkdir", line) and re.search('%{_exec_prefix}', line):
                print(re.sub("%{_exec_prefix}", "%{buildroot}/%{_exec_prefix}", line), file=fh)
            if re.search("^cp", line) and re.search('%{_exec_prefix}/bin', line):
                self.has_bin = True
                line = re.sub("^cp", "%{__install} -m 0755", line)
                line = re.sub("%{_exec_prefix}", "%{buildroot}/%{_exec_prefix}", line)
                print(line, file=fh)

    def parseprefix(self, line):
        line = re.sub("\\$PREFIX/bin", "%{_exec_prefix}/bin", line)
        line = re.sub("\\$PREFIX/lib", "%{_exec_prefix}/lib", line)
        line = re.sub("\\$PREFIX/include", "%{_prefix}/include", line)
        line = re.sub("\\$PREFIX", "%{_prefix}", line)
        return line

    def dump_file_content(self, fh):
        if self.has_bin:
            print("%{_exec_prefix}/bin/*", file=fh)

    def cleanup(self, fh):
        print("rm -fR %{_builddir}/%{name}-%{version}", file=fh)
        print("rm -fR %{buildroot}", file=fh)

    def download_source(self):
        # download the source file [[ count be moved to a function]]
        url = self.conda.conda_yaml["source"]["url"]
        a = urlparse(url)
        downto = os.path.join(self.args.rpm, "SOURCES", os.path.basename(a.path))

        if not os.path.exists(downto):
            logging.info(f"downloading")

            ctx = ssl.create_default_context()
            ctx.check_hostname = False
            ctx.verify_mode = ssl.CERT_NONE
            with urllib.request.urlopen(url, context=ctx) as u, open(downto, "wb") as file:
                file.write(u.read())

            logging.info(f"downloaded to [{downto}]")





conda2rpm_mapping = {"{{ compiler('cxx') }}": "gcc",
                     "{{ compiler('c') }}": "gcc",
                     "zlib": "zlib-devel", }


class SpecParser(Flounder):

    def __init__(self, spec):
        Flounder.__init__(self)
        self.spec = spec
        self.lines = None

    def get_spec_release(self):
        if self.lines is None:
            self.load_lines()
        version = None
        for line in self.lines:
            if re.search("^Release:", line):
                return re.search("(?<=Release: )[^%]+", line).group(0)
        return 1

    def load_lines(self):
        logging.info(f"reading SPEC file {self.spec}")
        with open(self.spec) as file:
            self.lines = file.readlines()


