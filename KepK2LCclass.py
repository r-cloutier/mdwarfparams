from imports import *
from truncate_cmap import *


def loadpickle(fname):
    fObj = open(fname, 'rb')
    self = pickle.load(fObj)
    fObj.close()
    return self


class KepK2LC:

    def __init__(self, name):
	try:
	    os.mkdir('PipelineResults')
	except OSError:
	    pass
	self.object_name = name
        try:
            os.mkdir('PipelineResults/%s'%self.object_name.replace(' ','_'))
        except OSError:
            pass
        self._pickleobject()


    def _pickleobject(self):
        fObj = open('%s/TESSResults'%self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close() 
