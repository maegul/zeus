import numpy as np

from circ_stat import circ_stat as cs

import Pickle
import glob
import copy


import collections

from . import hermes
from . import themis


class Athena:

	def __init__(self, name = None, description = None, path = '.'):

		'''
		path default is current directory.  Idea is all zeus files and athena files will be in same directory

		Decoupling of data and output path can be done at zeus level
		'''

		self.PROJ_ID = hermes.mk_project_ID(name = name, description = description)

		self.path = pthl.Path(path)

		# save absolute paths if necessary
		self._absolute_paths = {
			'path': self.path.absolute()
		}


	def save(self):

		file_name = self.path / ('Athena_' + self.PROJ_ID['name'] + '_' + self.PROJ_ID['ID'] + '.pkl')

		self.SavePath = file_name
		self._absolute_paths['SavePath'] = file_name.absolute()

		# Atomic saving (helpful?)
		temp_path = file_name.with_suffix(file_name.suffix + '.tmp')

		with open(temp_path, 'wb') as f:
			pickle.dump(self, f)

		temp_path.rename(file_name)

		print(f'Saved pickle file to {str(self.path)} as {str(file_name.name)}')

# def update(d, u):
# 	'''
# 	For updating nested dictionaries without destroying elements

# 	Parameters
# 	__________

# 	d: original dictionary
# 	u: dictionary to update the original by, it will be inserted where appropriate

# 	Returns
# 	_______
# 	d: original dictionary updated as is necessary

# 	'''

# 	for k, v in u.iteritems():

# 		# if branch is dict
# 		if isinstance(v, collections.Mapping):
# 			# go down branches of d that match u ...
# 			r = update(d.get(k, {}), v)
# 			d[k] = r
# 		else:
# 			d[k] = u[k]
# 	return d



# def meta_data_concat(spec, md):
#     '''
#     spec: list
#         specification of the meta data structure
#     md: list
#         actual meta data content to be used

#     Returns
#     _______
#     Concatenation of spec and md for use as keys in a nested dictionary
#     '''
#     return [param + '_' + str(md[param]) for param in spec]






# def metadata_to_dict(spec, md):
    
#     combination = meta_data_concat(spec, md)

#     d = reduce(lambda x, y: {y:x}, reversed(combination+[dict()]))

#     return d




# def add_data_leaf(proj, spec, md, key=None, obj=None):
    
#     dict_keys = meta_data_concat(spec, md)
    
#     dict_item = proj

#     for i, k in enumerate(dict_keys):
        
#         assert k in dict_item.keys(), '%s key is not in the dictionary' % k
        
#         dict_item = dict_item[k]
        
#         if i == len(dict_keys)-1:
#             dict_item.update({key: obj})
        
#     return proj




# def get_data_leaf(proj, spec, md, ext_keys=None):
    
#     dict_keys = meta_data_concat(spec, md)
    
#     if ext_keys is not None:
#         if type(ext_keys) == str:
#             ext_keys = [ext_keys]
            
#         dict_keys = dict_keys + ext_keys
    
#     dict_item = proj

#     for i, k in enumerate(dict_keys):
        
#         assert k in dict_item.keys(), '%s key is not in the dictionary' % k
        
#         dict_item = dict_item[k]
        
        
#     return dict_item




# def descs_type_test(descs):

#     assert type(descs) == list, 'Descriptions must be a list'
#     assert all([type(e) == str for e in descs]), 'Descriptions must be a list of strings'



# class athena:
#     '''
#     Defines the structure and contents of the meta-data for the project
    
#     Presumes a single analysis directory.
    
#     Saves a file for common reference to that directory.
    

    
#     '''

            
#     def __init__(self, project_name = None):
        
#         if project_name is not None:
#             self.project_name = project_name
            
#         else:

#         	# Can I just instantiate self as the pickled athena object
#         	# or iterate over each property and add to self?
#         	# really depends on what will be in the file?
#         		# prject meta data
#         		# all zeus instantiations
#         		# all tuning data and calculations
#         		# A nested ditionary?
#             panth_file_path = glob.glob('*.pantheon')
            
#             assert panth_file_path, 'There is no pantheon file for metadata upload'
            
#             assert len(panth_file_path) > 1, 'There are multiple pantheon files'
            
#             with open(panth_file_path[0], 'r') as f:
            
#                 self.proj_meta_data = cPickle.load(f)
                
#     def _set_proj(self, question, descs=None):
        
#         try:
#             self.proj_meta_data
            
#             print('meta data already loaded')
            
#         except AttributeError:
            
#             assert type(question) == str, 'question must be a string'
            
#             self.question = question

#             if descs is not None:
                
#                 descs_type_test(descs)
                
#                 descs = [e.lower() for e in descs]
                
#                 self.descs = descs
                
#     def _update_descs(self, update_descs, replace=False):
        
#         descs_type_test(update_descs)
        
#         try:
#             self.descs
            
#             update_descs = [e.lower() for e in update_descs]
            
#             if replace:
#                 self.descs = update_descs
#             else:
#                 for e in update_descs:
#                     self.descs.append(e)
                    
#         except AttributeError:
#             print('Must define descriptions before you can update')
            
            
    
#     def _zeus_meta_data(self, exp=None, track=None, unit=None, run=None, desc=None):
        
#         assert type(exp) == str, 'experiment should be a string'
#         assert type(track) == int, 'Track should be integer'
#         assert type(unit) == int, 'unit should be integer'
#         assert type(run) == int, 'run should be integer'
#         assert type(desc) == str, 'description should be a string'
        
#         exp = exp.upper()
#         desc = desc.lower()

#         self.__md_spec = ['experiment', 'track', 'unit', 'run']

        
#         meta_data_output = dict(project_name = self.project_name,
#                     question = self.question,
#                     experiment = exp.lower(),
#                     track = track,
#                     unit = unit,
#                     run = run
#                    )
    
#     	# Add description
#         try:
#             self.descs
            
#             if desc is not None:

#                 assert desc in self.descs, 'desc not prescribed'
                
#                 meta_data_output.update(desc=desc)

#         except AttributeError:
            
#             print('There are no prescribed descriptions ... make them!')

#         except AssertionError as e:
#             print('Your description is not prescribe')
#             print( repr(e))


#         try:
#             self.units = update(self.units, )


#         # use update to add unit metadata to master dictionary
#         # Add description to defaults?
#         # use add data to add pd cond_tuning to leaf
#         # Analysis functions on cond_tuning, using experiment type metadata from zeus? or just manual?

#         # use description to manage different experiments for each unit.
#         # cond_tuning stored in units.  Analysis done on that, and stored in unit.
#         # General, 'athena level' data is gathered from each experiment
#         # for each unit, all analysis results data is gathered to form a column in a dataframe
#         # each unit is a row, with a multi_index for the track and expeirment metadata also




            
#         meta_data_nested = dict(
#         	meta_data_output(self.__md_spec[0]) = dict(
#         		self.__md_spec[1] + str(meta_data_output(self.__md_spec[1])) = dict(
#         			self.__md_spec[2] + str(meta_data_output(self.__md_spec[2])) = dict(
#         				self.__md_spec[3] + str(meta_data_output(self.__md_spec[3])) = dict(
#         					desc = meta_data_output['desc']
#         					)
#         				)
#         			)
#         		)
#         	)


#         # Instantiate MetaData tree if it does not exist
#         try:
#         	assert type(self.units_meta_data) == dict, 'meta data structure should be a dictionary'

#     	except AttributeError:

# 	        self.units_meta_data = dict()




#         self.units_meta_data = update(self.units_meta_data, meta_data_nested)



            
#         return meta_data_output




#     def _add_data(self, zeus_instance):

#         try:
#             zeus_instance.cond_tuning

#             meta_data = copy.deepcopy(zeus_instance.meta_data)

#            	# attach whol zeus instance to nested dictionary?


            
#         except AttributeError:
            
#             print('zeus instance has no cond_tuning attribute ... analyse the data first!')
            

#         # Master list of all exp, tracks, units etc
#         # Check against this list for each new metadata ... checking if matches or not
#         # adding where new, with notification
#         # function for viewing the lists
            






# # Once metadata generated by athena, she remembers, and I can add with a simple funtion
# # The function strips the metadata from zeus and incorporates the (relevant) data
# # On the zeus side, all the unnecessary data is saved to file using a zeus method, with appropriate file name

 

