from time import time

def mkUID():
	'''
	Generates Unique ID from system time
	Returns hex of system epoch time
	'''

	return hex(int(time()*1e7))[2:]