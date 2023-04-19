"""
Created by: Fallon Ratner
Date: 12/04/23
Project: scRNAseq Prototype Analysis
"""
#import libraries
import pandas as pd
import numpy as np

#Read in output file from sctye_annot.R file

data1 = pd.read_csv("{Your file}",index_col=0)

data2 = pd.read_csv("{Your file}",index_col=0)

# Add a new 'Group' column to each DataFrame
data1['Group'] = '{Group Label}'
data2['Group'] = '{Group Label}'

#Calulate the proportion as a percentage
total_count = data1['ncells'].sum()

# Calculate the proportion of each cell type count
proportions = data1['ncells'] / total_count

# Add the proportions as a new column to the DataFrame
data1['Proportion'] = proportions

data1['Percentage'] = (proportions * 100)

sum = data1['Percentage'].sum()
# Display the updated DataFrame
print(data1)
#########
#Repeat for the 2nd df
total_count2 = data2['ncells'].sum()
# Calculate the proportion of each cell type count
proportions2 = data2['ncells'] / total_count2

# Add the proportions as a new column to the DataFrame
data2['Proportion'] = proportions2

data2['Percentage'] = (proportions2 * 100)
# Display the updated DataFrame
print(data2)

sum2 = data2['Percentage'].sum()

# Concatenate the two dataframes
combo = pd.concat([data2, data2])
combo.head()

#save new df as csv
combo.to_csv('{File name}.csv')

