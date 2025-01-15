import pandas as pd
from pymatgen.core import Composition
# Read the data
df = pd.read_csv("perovsk111_30kSamples_fullfeaturized_coverage.csv")

### This is the filter for exploration, convex hull coverage can be low for some compositions, 
### so one is advised to complement the exploration of the convex hull before calculating decomposition
### energies for some of these.
filtered_exploration_df = df[(df["stability"] <= 0.015) & (df["mlip_stability"] <= 0.045) &
                  (df["band_gap"] <= 3.5) & (df["gap_class"] >= 0.5) & (df["coverage_score"] >= 0.25)] 
### These structures are the ones considered in the diagrams for the overview of the method 
### presented in Fig.3 of the paper.

with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(filtered_exploration_df.to_string(index=False))

print(filtered_exploration_df.shape)
## Generate a csv
filtered_exploration_df.to_csv("filtered_exploration_df.csv",index=False)


## This is a more restricted filter which coverage_score is higher, therefore
## the decomposition energy for the structures can be estimated with more confidence.
filtered_production_df = df[(df["stability"] <= 0.015) & (df["mlip_stability"] <= 0.045) &
                  (df["band_gap"] <= 3.5) & (df["gap_class"] >= 0.5) & (df["coverage_score"] >= 0.45)]
## These are the structures that we consider for production with our workflow.

with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(filtered_production_df.to_string(index=False))
print(filtered_production_df.shape)
## Generate a csv
filtered_production_df.to_csv("filtered_production_df.csv",index=False)

def extract_Bele(structure):
    '''This is a function to organize the output by the B-site elements of the structure'''

    removeAlkaliMetals = ['Cs', 'Rb', 'K', 'Na', 'Li']
    removeHalogens = ['F', 'Cl', 'Br', 'I']
    removeChalcogens = ['O', 'S'] 
    chemicalSpace=Composition(structure).chemical_system.split('-')
    BchemicalSpace = [ele for ele in chemicalSpace if ele not in removeAlkaliMetals and 
                                                  ele not in removeHalogens and 
                                                  ele not in removeChalcogens]
    sortedBspace="-".join(sorted(BchemicalSpace))
    return sortedBspace

filtered_production_df['Bgroup'] = filtered_production_df['structure'].apply(lambda x: extract_Bele(x))
grouped_filtered_production_df = filtered_production_df.groupby('Bgroup') #.apply(lambda group: group.drop(columns=['initial_index'])).groupby('Bgroup')

def get_sorted_group(group,key):
     try:
          return group.get_group(key).sort_values(by='stability',ascending=True)
     except KeyError:
          print(f"No such key {key} in dictionary")

# just to exclude the initial index column that messes. 
with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(get_sorted_group(grouped_filtered_production_df,'Ag-Sb'))
    print(get_sorted_group(grouped_filtered_production_df,'Ag-Y'))
    print(get_sorted_group(grouped_filtered_production_df,'Co'))
