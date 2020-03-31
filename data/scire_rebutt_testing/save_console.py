
import shelve
import os
#T='Hiya'
#val=[1,2,3]

filename=os.getcwd()+'\\console_1'
my_shelf = shelve.open(filename,'n') # 'n' for new
variables = ['exp_dict_list_optimized','exp_dict_list_original']
for key in dir():
    if key in variables:
        print(key)
        my_shelf[key] = globals()[key]
    #except TypeError:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        #print('ERROR shelving: {0}'.format(key))
my_shelf.close()
#To restore:
