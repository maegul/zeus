from functools import wraps
import inspect

def getMethods(object_class, incl_hidden=False):

	methods = inspect.getmembers(
		object_class, predicate = inspect.isfunction
		)

	# method_strings = []
	# print('what?')

	# for m in methods:
	# 	# Method string -> ms
	# 	ms = m[0]

	# 	# If including hidden
	# 	if incl_hidden and ms[:2] == '__':
	# 		method_strings.append(ms)
	# 	elif ms[:2] != '__':
	# 		method_strings.append(ms)

	methods = [f[0] for f in methods]

	if not incl_hidden:
		# Use only not hidden methods
		methods = [f for f in methods if f[:2] != '__']


	# return method_strings
	
	return methods


def funcTracker(func):
	'''
	Track whether a function has been called or not

	Adds value to dictionary attribute of parent class
	'''

	@wraps(func)
	def wrapper(self, *args, **kwargs):

		func(self, *args, **kwargs)

		# presumes name of function_tracking stays the same
		# AND, that it has been initialised
		self._function_tracking.update(
				{
					func.__name__ : True
				}
			)

	return wrapper










