import json
import os
from collections import Counter
from functools import lru_cache
from importlib.resources import files
from itertools import combinations

import yaml

# Ancestors too generic to indicate a conflict: two classes that only meet here
# belong to different branches of the Biolink hierarchy.
BIOLINK_ROOTS = {'biolink:NamedThing', 'biolink:Entity'}

# Accepted class combinations, even when they share a parent
# (e.g. the Node Normalizer's gene/protein conflation).
ALLOWED_CLASS_COMBINATIONS = [
    {'biolink:Gene', 'biolink:Protein'},
]


def _to_curie(class_name):
    return 'biolink:' + ''.join(word[0].upper() + word[1:] for word in class_name.split())


@lru_cache(maxsize=None)
def load_biolink_hierarchy():
    """
    Load the Biolink model from the installed biolink-model package and return:
    - is_a: class -> direct parent (is_a only)
    - ancestors: class -> all ancestors (is_a + mixins)
    """
    schema_path = files('biolink_model').joinpath('schema/biolink_model.yaml')
    classes = yaml.safe_load(schema_path.read_text(encoding='utf-8'))['classes']

    is_a = {}
    parents = {}
    for name, definition in classes.items():
        definition = definition or {}
        curie = _to_curie(name)
        if definition.get('is_a'):
            is_a[curie] = _to_curie(definition['is_a'])
        parents[curie] = [_to_curie(p) for p in [definition.get('is_a')] + (definition.get('mixins') or []) if p]

    ancestors = {}
    def collect(curie):
        if curie not in ancestors:
            ancestors[curie] = set()
            for parent in parents.get(curie, []):
                ancestors[curie] |= {parent} | collect(parent)
        return ancestors[curie]
    for curie in parents:
        collect(curie)

    return is_a, ancestors


def _is_a_chain(curie, is_a):
    chain = [curie]
    while chain[-1] in is_a:
        chain.append(is_a[chain[-1]])
    return chain


@lru_cache(maxsize=None)
def find_category_conflicts(categories):
    """
    For a node's category set (frozenset), return:
    - the most specific classes (excluding classes that are ancestors of another)
    - conflicting pairs of specific classes and their nearest common is_a ancestor,
      when that ancestor is below NamedThing.
    """
    is_a, ancestors = load_biolink_hierarchy()
    most_specific = sorted(
        c for c in categories
        if not any(c in ancestors.get(other, ()) for other in categories if other != c)
    )

    conflicts = []
    for a, b in combinations(most_specific, 2):
        if {a, b} in ALLOWED_CLASS_COMBINATIONS:
            continue
        chain_b = _is_a_chain(b, is_a)
        common = next((x for x in _is_a_chain(a, is_a) if x in chain_b), None)
        if common is not None and common not in BIOLINK_ROOTS:
            conflicts.append({'classes': [a, b], 'common_ancestor': common})

    return most_specific, conflicts


def check_category_conflicts(input_file='data/nodes.jsonl', output_file='data/category_conflicts.jsonl'):
    """
    Flag nodes whose most specific classes descend from the same Biolink parent
    (a common is_a ancestor below NamedThing), e.g. MolecularMixture + SmallMolecule
    (common ancestor: ChemicalEntity). Classes from different branches are accepted.
    """
    line_count = 0
    flagged_count = 0
    flagged_combinations = Counter()

    try:
        with open(input_file, 'r', encoding='utf-8', newline='') as infile, \
             open(output_file, 'w', encoding='utf-8', newline='') as outfile:
            for line in infile:
                line_count += 1
                try:
                    node = json.loads(line)
                except json.JSONDecodeError:
                    print(f"\n[ATTENTION] Ligne {line_count} mal formatée (skip).")
                    continue

                categories = node.get('category', [])
                if isinstance(categories, str):
                    categories = [categories]
                most_specific, conflicts = find_category_conflicts(frozenset(categories))

                if conflicts:
                    flagged_count += 1
                    flagged_combinations[tuple(most_specific)] += 1
                    # Preserve the node as-is by writing its original input line.
                    outfile.write(line if line.endswith('\n') else line + '\n')
    except FileNotFoundError:
        print(f"FATAL ERROR : Le fichier d'entrée '{input_file}' n'a pas été trouvé.")
        return

    print(f"\n--- Category Conflict Check Complete ---")
    print(f"Total lines read: {line_count}")
    print(f"Nodes with conflicting most specific categories: {flagged_count}")
    for combination, count in flagged_combinations.most_common():
        print(f"  {count:>8}  {' + '.join(c.removeprefix('biolink:') for c in combination)}")
    print(f"Successfully exported results to: {output_file}")


def filter_and_export_nodes(input_file='data/nodes.jsonl', output_file='problematic_biolink_nodes.jsonl',REQUIRED_CATEGORIES = {'biolink:MolecularMixture', 'biolink:SmallMolecule'}):
    """
    Read nodes.jsonl, find nodes that have ALL required categories
    ('biolink:MolecularMixture' AND 'biolink:SmallMolecule'), and export them.
    """
    
    # Both required categories (AND logic)
    catchecks = len(REQUIRED_CATEGORIES)
    
    filtered_nodes = []
    filtered_nodes_intersection = []
    no_category_nodes = []
    
    print(f"--- Démarrage du Filtrage ---")

    try:
        with open(input_file, 'r', encoding='utf-8') as infile:
            line_count = 0
            found_count = 0
            found_count_intersection = 0
            found_count_no_category = 0
            
            for line in infile:
                line_count += 1
                try:
                    node = json.loads(line)
                    categories_raw = node.get('category',[])

                    if len(categories_raw) == 0:
                        print(f"\n[ATTENTION] Ligne {line_count} : no category.")
                        no_category_nodes.append(node)
                        found_count_no_category += 1

                    current_node_categories = set()
                    if isinstance(categories_raw, str): # if len ==1
                        current_node_categories.add(categories_raw)
                    elif isinstance(categories_raw, list):
                        for cat in categories_raw:
                            if isinstance(cat, str):
                                current_node_categories.add(cat)
                            elif isinstance(categories_raw, list):
                                error_msg = f"\n[ATTENTION] Ligne {line_count} : catégorie mal formatée (skip)."
                                print(error_msg)

                    if REQUIRED_CATEGORIES.issubset(current_node_categories):
                        filtered_nodes.append(node)
                        found_count += 1

                except json.JSONDecodeError:
                    print(f"\n[ATTENTION] Ligne {line_count} mal formatée (skip).")
                    continue

            for filter_node in filtered_nodes:
                if len(filter_node.get('category', [])) == catchecks:
                    filtered_nodes_intersection.append(filter_node)
                    found_count_intersection += 1


    except FileNotFoundError:
        print(f"FATAL ERROR : Le fichier d'entrée '{input_file}' n'a pas été trouvé.")
        return

    # Export:
    try:
        with open(os.path.join(os.getcwd(), output_file), 'w', encoding='utf-8') as outfile:
            for node in filtered_nodes:
                outfile.write(json.dumps(node) + '\n')
    except IOError as e:
        print(f"\nFATAL ERROR lors de l'écriture du fichier résultat : {e}")
        return

    # Final Summary
    total_lines = line_count
    if total_lines > 0:
        percentage = (found_count / total_lines) * 100
        print(f"\n--- Processing Complete ---")
        print(f"Total lines read: {total_lines}")
        print(f"Target nodes found: {found_count}")
        print(f"Filter rate: {percentage:.2f}%")
        print(f"Nodes with exactly {catchecks} categories: {found_count_intersection}")
        print(f"Nodes with no categories: {found_count_no_category}")
        print(f"Nodes with at least {catchecks} of the required categories: {output_file}")
        print(f"Successfully exported results to: {output_file}")

# Run the function
if __name__ == "__main__":
    filter_and_export_nodes()
