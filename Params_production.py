import numpy as np
import json 
import os
import pandas as pd
import datetime
import sys

from params_values import Dictionary as Dictionary

Max_Time=(pd.to_datetime(Dictionary['Last_colocation_file'])-pd.to_datetime(Dictionary['time_starting_observation'])).days+7
Number_of_colocation=int(Max_Time/7)
with open(sys.argv[3]+'/dict_polygon_name.txt') as f: 
    dict_polygon_name = f.read()
dict_polygon_name = json.loads(dict_polygon_name) 
dict_polygon_name=dict((key, value-1) for (key, value) in dict_polygon_name.items())
Seeding_province=dict_polygon_name[sys.argv[1]]
Instance_of_infection=(pd.to_datetime(sys.argv[2])-pd.to_datetime(Dictionary['time_starting_observation'])).days
float_numbers=[Dictionary["alpha"],Dictionary["mu"]]
np.savetxt("params_double_"+sys.argv[1]+".txt", float_numbers)
integer_numbers=[Dictionary["Tiles"],Max_Time,Number_of_colocation,Instance_of_infection,Seeding_province,Dictionary["Initial_infected"]]
np.savetxt("params_int_"+sys.argv[1]+".txt", integer_numbers, fmt="%d")
