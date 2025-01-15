import mysql.connector
import pandas as pd
from pymatgen.core import Composition
from itertools import combinations
import os
import logging
import json

# --- Logging Configuration ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- Database Configuration ---
DB_CONFIG = {
    'user': 'rogerio',
    'password': '123456',
    'host': '%',
    'port': 3306,
    'db': 'oqmd_db',
    'unix_socket': '/var/run/mysqld/mysqld.sock'
}

# --- Lookup Table for Total Combinations ---
LOOKUP_TABLE = {
    2: {'n_pairs': 1, 'n_triples': 0, 'n_quadruples': 0, 'n_quintuples': 0, 'n_sextuples': 0, 'n_septuples': 0, 'n_octuples': 0},
    3: {'n_pairs': 3, 'n_triples': 1, 'n_quadruples': 0, 'n_quintuples': 0, 'n_sextuples': 0, 'n_septuples': 0, 'n_octuples': 0},
    4: {'n_pairs': 6, 'n_triples': 4, 'n_quadruples': 1, 'n_quintuples': 0, 'n_sextuples': 0, 'n_septuples': 0, 'n_octuples': 0},
    5: {'n_pairs': 10, 'n_triples': 10, 'n_quadruples': 5, 'n_quintuples': 1, 'n_sextuples': 0, 'n_septuples': 0, 'n_octuples': 0},
    6: {'n_pairs': 15, 'n_triples': 20, 'n_quadruples': 15, 'n_quintuples': 6, 'n_sextuples': 1, 'n_septuples': 0, 'n_octuples': 0},
    7: {'n_pairs': 21, 'n_triples': 35, 'n_quadruples': 35, 'n_quintuples': 21, 'n_sextuples': 7, 'n_septuples': 1, 'n_octuples': 0},
    8: {'n_pairs': 28, 'n_triples': 56, 'n_quadruples': 70, 'n_quintuples': 56, 'n_sextuples': 28, 'n_septuples': 8, 'n_octuples': 1}
}

def fetch_data_from_database():
    """Fetches data from the MySQL database."""
    logging.info("Fetching data from the MySQL database...")
    cnx = mysql.connector.connect(**DB_CONFIG)
    cursor = cnx.cursor()
    query = "SELECT id, composition_id, entry_id, stability FROM formation_energies WHERE stability < 2.89"
    cursor.execute(query)
    data = cursor.fetchall()
    cnx.close()
    logging.info("Data fetched from the database.")
    return pd.DataFrame(data, columns=['id', 'composition_id', 'entry_id', 'stability'])

def remove_polymorphs(df):
    """Removes polymorphs by keeping the first entry for each composition."""
    logging.info("Removing polymorphs...")
    df = df.sort_values('stability').groupby('composition_id').first().reset_index()
    logging.info("Polymorphs removed.")
    return df

def normalize_composition(composition_str):
    """Normalizes a composition string by removing stoichiometry."""
    comp = Composition(composition_str)
    return comp.reduced_formula

def reduce_composition_to_elements(composition_str):
    """Transforms a composition string into individual elements with '1' stoichiometry."""
    comp = Composition(composition_str)
    return " ".join([f"{el}1" for el in comp.elements])

def create_one_hot_encoded_df(df):
    """Creates a DataFrame with one-hot encoded elements."""
    logging.info("Creating one-hot encoded DataFrame...")
    df['normalized_composition'] = df['composition_id'].apply(normalize_composition)
    df['reduced_composition'] = df['normalized_composition'].apply(reduce_composition_to_elements)
    df = df.drop_duplicates(subset=['reduced_composition'])
    
    compositions = [Composition(comp) for comp in df['reduced_composition']]
    all_elements = set(el for composition in compositions for el in composition.elements)
    
    one_hot_encoded_data = []
    for composition in compositions:
        one_hot_vector = {
            'composition': str(composition),
            **{str(el): (el in composition.elements) for el in all_elements}
        }
        one_hot_encoded_data.append(one_hot_vector)
    
    one_hot_df = pd.DataFrame(one_hot_encoded_data).fillna(0)
    for col in one_hot_df.columns:
        if col != 'composition':
            one_hot_df[col] = one_hot_df[col].astype(int)
    logging.info("One-hot encoded DataFrame created.")
    return one_hot_df

def analyze_composition_coverage(composition_str, df):
    """
    Analyze the coverage of element combinations in the database.
    Only counts first exact match for each combination and accounts for total possible combinations.
    """
    comp = Composition(composition_str)
    elements = list(comp.elements)
    n_elements = len(elements)
    
    # Get total combinations from the lookup table
    total_combinations = LOOKUP_TABLE.get(n_elements, {
        'n_pairs': 0, 'n_triples': 0, 'n_quadruples': 0, 
        'n_quintuples': 0, 'n_sextuples': 0, 'n_septuples': 0, 'n_octuples': 0
    })
    
    # Initialize counts
    counts = {
        'n_pairs': 0, 'n_triples': 0, 'n_quadruples': 0,
        'n_quintuples': 0, 'n_sextuples': 0, 'n_septuples': 0, 'n_octuples': 0,
        'total_combinations': total_combinations
    }
    
    # Get all elements in the database (excluding 'composition' column)
    all_elements = [col for col in df.columns if col != 'composition']
    
    if n_elements >= 2:
        for pair in combinations(elements, 2):
            pair = [str(el) for el in pair]
            pair_mask = df[pair].all(axis=1)
            other_elements = [el for el in all_elements if el not in pair]
            no_other_elements = ~df[other_elements].any(axis=1)
            exact_pair_mask = pair_mask & no_other_elements
            if exact_pair_mask.any():
                counts['n_pairs'] += 1
    
    if n_elements >= 3:
        for triple in combinations(elements, 3):
            triple = [str(el) for el in triple]
            triple_mask = df[triple].all(axis=1)
            other_elements = [el for el in all_elements if el not in triple]
            no_other_elements = ~df[other_elements].any(axis=1)
            exact_triple_mask = triple_mask & no_other_elements
            if exact_triple_mask.any():
                counts['n_triples'] += 1
    
    if n_elements >= 4:
        for quad in combinations(elements, 4):
            quad = [str(el) for el in quad]
            quad_mask = df[quad].all(axis=1)
            other_elements = [el for el in all_elements if el not in quad]
            no_other_elements = ~df[other_elements].any(axis=1)
            exact_quad_mask = quad_mask & no_other_elements
            if exact_quad_mask.any():
                counts['n_quadruples'] += 1
    
    if n_elements >= 5:
        for quintuple in combinations(elements, 5):
            quintuple = [str(el) for el in quintuple]
            quintuple_mask = df[quintuple].all(axis=1)
            other_elements = [el for el in all_elements if el not in quintuple]
            no_other_elements = ~df[other_elements].any(axis=1)
            exact_quintuple_mask = quintuple_mask & no_other_elements
            if exact_quintuple_mask.any():
                counts['n_quintuples'] += 1
    
    if n_elements >= 6:
        for sextuple in combinations(elements, 6):
            sextuple = [str(el) for el in sextuple]
            sextuple_mask = df[sextuple].all(axis=1)
            other_elements = [el for el in all_elements if el not in sextuple]
            no_other_elements = ~df[other_elements].any(axis=1)
            exact_sextuple_mask = sextuple_mask & no_other_elements
            if exact_sextuple_mask.any():
                counts['n_sextuples'] += 1
    
    if n_elements >= 7:
        for septuple in combinations(elements, 7):
            septuple = [str(el) for el in septuple]
            septuple_mask = df[septuple].all(axis=1)
            other_elements = [el for el in all_elements if el not in septuple]
            no_other_elements = ~df[other_elements].any(axis=1)
            exact_septuple_mask = septuple_mask & no_other_elements
            if exact_septuple_mask.any():
                counts['n_septuples'] += 1
    
    if n_elements >= 8:
        for octuple in combinations(elements, 8):
            octuple = [str(el) for el in octuple]
            octuple_mask = df[octuple].all(axis=1)
            other_elements = [el for el in all_elements if el not in octuple]
            no_other_elements = ~df[other_elements].any(axis=1)
            exact_octuple_mask = octuple_mask & no_other_elements
            if exact_octuple_mask.any():
                counts['n_octuples'] += 1
    
    # Calculate coverage percentages
    counts['coverage'] = {
        'pairs': (counts['n_pairs'] / total_combinations['n_pairs']) if total_combinations['n_pairs'] > 0 else 0,
        'triples': (counts['n_triples'] / total_combinations['n_triples']) if total_combinations['n_triples'] > 0 else 0,
        'quadruples': (counts['n_quadruples'] / total_combinations['n_quadruples']) if total_combinations['n_quadruples'] > 0 else 0,
        'quintuples': (counts['n_quintuples'] / total_combinations['n_quintuples']) if total_combinations['n_quintuples'] > 0 else 0,
        'sextuples': (counts['n_sextuples'] / total_combinations['n_sextuples']) if total_combinations['n_sextuples'] > 0 else 0,
        'septuples': (counts['n_septuples'] / total_combinations['n_septuples']) if total_combinations['n_septuples'] > 0 else 0,
        'octuples': (counts['n_octuples'] / total_combinations['n_octuples']) if total_combinations['n_octuples'] > 0 else 0
    }
    
    return counts

def calculate_coverage_score(coverage):
    """
    Calculate a weighted score based on the coverage percentages, only considering combinations
    up to the number of elements in the composition.
    Weights are assigned as follows: 1 for pairs, 3 for triples, 6 for quadruples, 10 for quintuples,
    15 for sextuples, 21 for septuples, and 28 for octuples.
    The different weights are there to reflect that if you have information
    in a structure that shares more elements to the target composition, it is more valuable than
    lower level combinations.
    """
    weights_per_class = {
        'pairs': 1,
        'triples': 3,
        'quadruples': 6,
        'quintuples': 10,
        'sextuples': 15,
        'septuples': 21,
        'octuples': 28
    }
    
    weighted_score = 0
    total_weight = 0
    
    # Map total_combinations keys to coverage keys
    mapping = {
        'n_pairs': 'pairs',
        'n_triples': 'triples',
        'n_quadruples': 'quadruples',
        'n_quintuples': 'quintuples',
        'n_sextuples': 'sextuples',
        'n_septuples': 'septuples',
        'n_octuples': 'octuples'
    }
    
    # Only include weights for combinations that exist in the composition
    for total_key, coverage_key in mapping.items():
        if coverage['total_combinations'][total_key] > 0:
            weighted_score += weights_per_class[coverage_key] * coverage['coverage'][coverage_key]
            total_weight += weights_per_class[coverage_key]
    
    return weighted_score / total_weight if total_weight > 0 else 0

def main(csv_file='CSIp3m1_30kselected_relaxed_featurized_OmegaROSA_predictions.csv', 
         oqmd_csv_file='OQMD_one_hot_encoded_data.csv', 
         output_csv='CSIp3m1_30kselected_relaxed_featurized_OmegaROSA_predictions_coverage.csv',
         progress_file='progress.json',
         chunk_size=10):
    
    logging.info(f"Starting the incremental process with chunk size: {chunk_size}")
    
    # Load the selected.csv file
    selected_df = pd.read_csv(csv_file)
    logging.info(f"Loaded selected DataFrame from: {csv_file}")
    
    # Load or create one-hot encoded dataframe
    if os.path.exists(oqmd_csv_file):
        one_hot_df = pd.read_csv(oqmd_csv_file)
        logging.info(f"Loaded existing one-hot encoded DataFrame from: {oqmd_csv_file}")
    else:
        # Fetch data from the database and process it
        df = fetch_data_from_database()
        df = remove_polymorphs(df)
        
        # Create one-hot encoded dataframe
        one_hot_df = create_one_hot_encoded_df(df)
        logging.info(f"Created new one-hot encoded DataFrame")
        
        # Save the one-hot encoded DataFrame to a CSV file
        one_hot_df.to_csv(oqmd_csv_file, index=False)
        logging.info(f"One-hot encoded DataFrame saved to: {oqmd_csv_file}")
    
    # Load progress from file, if it exists
    start_index = 0
    if os.path.exists(progress_file):
        with open(progress_file, 'r') as f:
            progress_data = json.load(f)
            start_index = progress_data.get('last_processed_index', 0)
        logging.info(f"Resuming from index: {start_index}")
    else:
        logging.info("Starting from the beginning.")
    
    # Process in chunks
    for i in range(start_index, len(selected_df), chunk_size):
        end_index = min(i + chunk_size, len(selected_df))
        chunk_df = selected_df.iloc[i:end_index]
        logging.info(f"Processing chunk from index {i} to {end_index}")
        
        coverage_scores = []
        for j, composition in enumerate(chunk_df['structure']):
            logging.info(f"Processing composition {j+1}/{len(chunk_df)}: {composition}")
            coverage = analyze_composition_coverage(composition, one_hot_df)
            coverage_score = calculate_coverage_score(coverage)
            coverage_scores.append(coverage_score)
            logging.info(f"Coverage score calculated: {coverage_score}")
        
        chunk_df['coverage_score'] = coverage_scores
        
        # Append to output CSV
        if os.path.exists(output_csv):
            chunk_df.to_csv(output_csv, mode='a', header=False, index=False)
        else:
            chunk_df.to_csv(output_csv, index=False)
        
        # Save progress
        with open(progress_file, 'w') as f:
            json.dump({'last_processed_index': end_index}, f)
        logging.info(f"Progress saved. Last processed index: {end_index}")
    
    logging.info(f"Incremental process completed. Results saved to: {output_csv}")

if __name__ == "__main__":
    main()
