#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert pickle file to csv for transferal to R.

@author: murph
"""

import pandas as pd
import pickle
import os


os.chdir('/Users/murph/GitLab/sa_for_chemistry/fixed_dfn/data')

data = pd.read_pickle('single_fracture_combined_data.pkl')

# with open("combined_data.pkl", "rb") as pfile:
#     data = pickle.load(pfile) 

data.to_csv("combined_data.csv")
