from imports import *
from truncate_cmap import *


def loadpickle(fname):
    fObj = open(fname, 'rb')
    try:
        self = pickle.load(fObj)
    except UnicodeDecodeError:
        # necessary to load python v2.7 pickles in v3.6 
        self = pickle.load(fObj, encoding='latin1')
    fObj.close()
    return self


class LCclass:

    def __init__(self, folder, name, index):
        try:
            os.mkdir(folder)
        except OSError:
            pass
        self.object_name = str(name)
        self.folder_full = '%s/%s'%(folder, self.object_name.replace(' ','_'))
        self.fname_full = '%s/LC_%.5d'%(self.folder_full, index)
        try:
            os.mkdir(self.folder_full)
        except OSError:
            pass
        self._pickleobject()


    def _pickleobject(self):
        fObj = open(self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close() 
