import category_checks


def main():
    print("Hello from kg-biolink-sanity-check!")
    category_checks.check_category_conflicts(input_file='data/nodes.jsonl', output_file='data/category_conflicts.jsonl')


if __name__ == "__main__":
    main()
