import cPickle as pickle


class catalog_object:


    def __init__(self, fname):
	self.fname = fname


    def _pickleobject(self):
	f = open(self.fname, 'wb')
	pickle.dump(self, f)
	f.close()


def loadpickle(fname):
    f = open(fname, 'rb')
    self = pickle.load(f)
    f.close()
    return self

