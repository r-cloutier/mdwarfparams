from imports import *
from truncate_cmap import *


def loadpickle(fname):
    fObj = open(fname, 'rb')
    self = pickle.load(fObj)
    fObj.close()
    return self


class TESSLC:

    def __init__(self, TICid, index):
	try:
	    os.mkdir('PipelineResults')
	except OSError:
	    pass
	self.TICid = TICid
	self.folder_full = 'PipelineResults/TIC_%i'%self.TICid
        self.fname_full = '%s/TESSLC_%.5d'%(self.folder_full, index)
	try:
            os.mkdir(self.folder_full)
        except OSError:
            pass
        self._pickleobject()


    def _pickleobject(self):
        fObj = open(self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close() 
