import gzip
import os
import re
import tqdm


class scan_content:
    def __init__(self, provided_path, regex):
        print("searching [%s] for target files" % provided_path)
        self.regex = regex
        self.candidates = []
        self.scan(provided_path)
        self.report()

    def scan(self, path):
        files = os.listdir(path)
        for f in files:
            qfile = os.path.join(path, f)
            if os.path.isdir(qfile):
                self.scan(qfile)
            else:
                if re.search(self.regex, f):
                    self.define_hit(qfile)

    def define_hit(self, cpath):
        # print("identified candidate file [%s]" % cpath)
        self.candidates.append(cpath)

    def report(self):
        print("%s candidate file(s) found" % self.candidate_count())

    def candidate_count(self):
        return len(self.candidates)

    def aggregate(self, target):
        with gzip.GzipFile(target, 'wb') as fout:
            for f in tqdm(self.candidates):
                # print(f)
                file = open(f, "r")
                fout.write(file.read().encode('utf-8'))
                file.close()
            fout.close()
