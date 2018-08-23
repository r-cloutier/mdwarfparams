from imports import *

class GAIAMdwarfs():

    def __init__(self, fname):
	self.fname = fname


    def _pickleoject(self, fname=''):
	fname = self.fname if fname == '' else fname
        fObj = open('%s'%fname, 'wb')
        pickle.dump(self, fObj)
        fObj.close()
