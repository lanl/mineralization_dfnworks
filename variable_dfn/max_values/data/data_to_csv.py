#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert pickle file to csv for transferal to R.

@author: murph
"""

import pandas as pd
import pickle
import os


os.chdir('/Users/murph/GitLab/sa_for_chemistry/variable_dfn/max_values/data')

data = pd.read_pickle('multi_fracture_combined_data_max.pkl')

# with open("combined_data.pkl", "rb") as pfile:
#     data = pickle.load(pfile) 

data.to_csv("combined_data.csv")
