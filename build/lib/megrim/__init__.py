import os
import progressbar
import unicodedata
import urllib

from .extract_file_system_content import extract_file_system_content
from .project_manager import project_manager
from .project_registry import project_registry
from .project_selector import project_selector
from .reference_data_registry import reference_data_registry
from .scan_content import scan_content


def web_check(url):
    try:
        thepage = urllib.request.urlopen(url)
    except HTTPError as e:
        return False
    except URLError as e:
        return False
    else:
        return True


def removeDisallowedFilenameChars(filename):
    # https://stackoverflow.com/questions/295135/turn-a-string-into-a-valid-filename
    validFilenameChars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    cleanedFilename = unicodedata.normalize(u'NFKD', filename).encode('ASCII', 'ignore').decode('ascii')
    return ''.join(c for c in cleanedFilename if c in validFilenameChars)


def target_dir(target):
    if not os.path.exists(target):
        try:
            os.makedirs(target)
            print('directory [%s] created' % target)
            return True
        except:
            print('unable to create target [%s]' % target)
            raise RuntimeError('unable to create target [%s]' % xtarget)
    return False


class MyProgressBar():
    def __init__(self):
        self.pbar = None

    def __call__(self, block_num, block_size, total_size):
        if not self.pbar:
            self.pbar = progressbar.ProgressBar(maxval=total_size)
            self.pbar.start()

        downloaded = block_num * block_size
        if downloaded < total_size:
            self.pbar.update(downloaded)
        else:
            self.pbar.finish()
