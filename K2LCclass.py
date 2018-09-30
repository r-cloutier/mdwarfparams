from imports import *
from truncate_cmap import *


def loadpickle(fname):
    fObj = open(fname, 'rb')
    self = pickle.load(fObj)
    fObj.close()
    return self


class K2LC:

    def __init__(self, name, index):
	try:
	    os.mkdir('PipelineResults')
	except OSError:
	    pass
	self.object_name = name
	self.folder_full = 'PipelineResults/%s'%self.object_name.replace(' ','_')
        self.fname_full = '%s/K2LC_%.5d'%(self.folder_full, index)
	try:
            os.mkdir(self.folder_full)
        except OSError:
            pass
        self._pickleobject()


    def _pickleobject(self):
        fObj = open(self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close() 
