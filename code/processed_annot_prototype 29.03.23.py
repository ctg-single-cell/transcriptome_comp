"""
Created by: Fallon Ratner
Date: 29/03/23
Project: scRNAseq Prototype Analysis
"""

#import libraries
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

#file path: C:\Users\fallo\OneDrive\Documents\Internship_2023\transcriptome_comp\output


#read in Velsco files: processed data
vel = pd.read_csv('counts_3and6mon.csv')

#make datframe for 3months only 
mon3 = vel[vel['Group'].str.contains('3months')]

#make datframe for 6months only
mon6 = vel[vel['Group'].str.contains('6months')]
#############################################################################
#Make a column for proportion in mon3 dataframes
# Calculate the sum of all cell counts
total_count = mon3['Count'].sum()

# Calculate the proportion of each cell type count
proportions = mon3['Count'] / total_count

# Add the proportions as a new column to the DataFrame
mon3['Proportion'] = proportions

# Convert the proportion column to a percentage
mon3['Percentage'] = (proportions * 100).apply(np.ceil).astype(int)

# Display the updated DataFrame
print(mon3)

#Remove the first column
mon3.pop(mon3.columns[0])

#Save new df as csv
mon3.to_csv('proportion_velasco_3mon.csv')

#############################################################################
#Make a column for proportion in mon6 dataframes
# Calculate the sum of all cell counts
total_counts = mon6['Count'].sum()

# Calculate the proportion of each cell type count
proportions2 = mon6['Count'] / total_counts

# Add the proportions as a new column to the DataFrame
mon6['Proportion'] = proportions2

# Convert the proportion column to a percentage
mon6['Percentage'] = (proportions2 * 100).apply(np.ceil).astype(int)

# Display the updated DataFrame
print(mon6)

#Remove the first column
mon6.pop(mon6.columns[0])

#Save new df as csv
mon6.to_csv('proportion_velasco_6mon.csv')

#################################################################################
#Read in combined Velsco & Polioudakis file: processed data
combo = pd.read_csv('proportion_velasco_polio.csv')

#remove the first column
combo.pop(combo.columns[0])

#################################################################################
#Plot only IP (intermediate progenitors)
ip = combo[combo['CellType'].str.contains('IP')]

test =pd.DataFrame(ip, columns=('CellType', 'Group', 'Percentage'))

combo.to_csv('proportion_velasco_poliov2.csv')

# Create the stacked bar chart
ax = test.plot(kind='bar')

# Add labels and title 
ax.set_xticklabels(labels, fontsize = 7)
ax.set_xlabel('Cell Type')
ax.set_ylabel('Proportion (%)')
ax.set_title('Intermediate Progenitors Percentage by Group Velasco et al., 2019 & Polioudakis et al., 2019')


# Create the stacked bar chart with seaborn
ax = sns.barplot(x='CellType', y='Percentage', data=test, hue='Group', dodge=False)
plt.title("Intermediate Progenitors Percentage by Group Velasco et al., 2019 & Polioudakis et al., 2019")
plt.xticks(rotation=90)
plt.ylabel("Proportion (%)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

##############################################################################
#Plot only Cycling/Pg2M/PgS
#In Polioudakis combine Pg2M and PgS values & rename it Cycling
cy = combo[combo['CellType'].str.contains('Pg')]
#sum cycling variables
new_row = cy.sum(axis=0)
combo2 = combo.append(new_row, ignore_index=True)
combo2['CellType'] = combo2['CellType'].replace('PgSPgG2M', 'Cycling')

cy2 = combo2[combo2['CellType'].str.contains('Cycling')]

test2 =pd.DataFrame(cy2, columns=('CellType', 'Group', 'Percentage'))
test2['Group'] = test2['Group'].replace('17-18 GW17-18 GW', '17-18 GW')

#Save df as csv
test2.to_csv('proportion_velasco_polio_cycling.csv')

# Create the stacked bar chart
ax = test2.plot(kind='bar')

# Add labels and title
ax.set_xlabel('Cell Type')
ax.set_ylabel('Proportion (%)')
ax.set_title('Cycling Cells: Velasco et al., 2019 & Polioudakis et al., 2019')


# Create the stacked bar chart with seaborn
ax = sns.barplot(x='CellType', y='Percentage', data=test2, hue='Group')
plt.title("Cycling Cells: Velasco et al., 2019 & Polioudakis et al., 2019")
plt.xticks(rotation=90)
plt.ylabel("Proportion (%)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
