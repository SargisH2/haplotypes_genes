import yaml
import re

with open('data/genes.yaml', 'r') as file:
    original_data = yaml.safe_load(file)
    
def get_substring_coordinates(main_string: str, substring: str) -> list:
    """
    Finds all occurrences of a substring within a string and returns their start and end indices.

    Parameters:
    - main_string (str): The string in which to search for the substring.
    - substring (str): The substring to find within the main_string.

    Returns:
    - list: A list of lists, where each sublist contains the start and end indices (inclusive, exclusive) of each occurrence of the substring.
    """
    coordinates = []
    start = 0
    
    main_string = main_string.lower()
    substring = substring.lower()
    while True:
        start = main_string.find(substring, start)
        if start == -1:
            break
        end = start + len(substring)
        coordinates.append([start, end])
        start += 1
    
    return coordinates

def check_gene(input_str: str, gene: dict) -> dict | None:
    """
    Checks if a gene is present in the input string and returns its information if found.

    Parameters:
    - input_str (str): The input string to search for the gene.
    - gene (dict): A dictionary containing information about the gene, including its synonyms.

    Returns:
    - dict or None: A dictionary containing the gene information if found, or None if not found.
    """
    all_gene_coords = []
    synonyms_to_check = gene.get('synonyms', [])
    if 0 == sum([gene['name'] in synonym for synonym in synonyms_to_check]):
        synonyms_to_check += [gene['name']]
    for synonym in synonyms_to_check:
        synonym_coords = get_substring_coordinates(input_str, synonym)
        for synonym_coord in synonym_coords: # [[substr1_coords], [substr1_other_coords]]
            if synonym_coord[0] not in [coord[0] for coord in all_gene_coords]:
                all_gene_coords.extend(synonym_coords)
            else:
                for coord in all_gene_coords:
                    if coord[0] == synonym_coord[0]:
                        coord[1] = max(coord[1], synonym_coord[1])
    if all_gene_coords:
        gene_found = {
            'name': gene.get('name'),
            'positions': all_gene_coords
        }
        return gene_found
        
def get_genes(input_str: str) -> dict:
    """
    Finds genes present in the input string and returns their information.

    Parameters:
    - input_str (str): The input string to search for genes.

    Returns:
    - dict: A dictionary containing information about genes found in the input string.
    """
    genes_found = []
    for gene in original_data:
        result = check_gene(input_str, gene)
        if result is not None:
            genes_found.append(result)
    return {
        "genes": genes_found
    }
    

def get_hla_match(hla_string: str) -> dict | None:
    """
    Parses an HLA string to extract components based on several regex patterns.

    Parameters:
    - hla_string (str): The string representing an HLA type.

    Returns:
    - dict: A dictionary containing the 'gene', 'allele_group', and 'protein' components of the HLA string.
      Returns None for 'allele_group' and 'protein' if they are not applicable.
      Returns None if no matches
    """
    # cases with concatenated allele group and protein
    pattern = r"HLA-([A-Z]+[0-9]*)(?:\*|:|\.)(\d{2})(\d{2})?"
    # without delimiters between allele and protein
    pattern_compact = r"HLA-([A-Z]+)(\d)\.(\d)"
    # delimiters between allele and protein
    pattern_delimited = r"HLA-([A-Z]+[0-9]*)(?:\*|:|\.)(\d{2})(?:\*|:|\.)(\d{2})"
    # gene only
    pattern_gene_only = r"HLA-([A-Z]+[0-9]*)$"

    match_compact = re.search(pattern_compact, hla_string)
    if match_compact:
        return {
            "gene": match_compact.group(1),
            "allele": match_compact.group(2),
            "protein": match_compact.group(3)
        }

    match_delimited = re.search(pattern_delimited, hla_string)
    if match_delimited:
        return {
            "gene": match_delimited.group(1),
            "allele": match_delimited.group(2),
            "protein": match_delimited.group(3)
        }

    match = re.search(pattern, hla_string)
    if match:
        return {
            "gene": match.group(1),
            "allele": match.group(2),
            "protein": match.group(3) if match.group(3) else None
        }

    match_gene_only = re.search(pattern_gene_only, hla_string)
    if match_gene_only:
        return {
            "gene": match_gene_only.group(1),
            "allele": None,
            "protein": None
        }

def get_hla(input_str: str) -> dict:
    """
    Identifies HLA markers in a provided string, extracting them along with their position within the string.

    Parameters:
    - input_str (str): The string from which to identify and extract HLA markers.

    Returns:
    - dict: A dictionary with a key 'hla' containing a list of dictionaries for each identified HLA marker.
      Each dictionary includes the HLA marker components and the positions of the marker within the input string.
    """
    hla_found = []
    text_words = [word.strip('.,') for word in input_str.split()]
    word_match_tuples = []
    for word in text_words:
        word_match = get_hla_match(word)
        if word_match: word_match_tuples.append((word, word_match))
    
    for hla_raw, match_dict in word_match_tuples:
        match_dict.update({
            "positions": get_substring_coordinates(input_str, hla_raw)
        })
        hla_found.append(match_dict)
    
    return {
        'hla': hla_found
    }