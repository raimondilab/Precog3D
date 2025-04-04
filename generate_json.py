import argparse
import json
import os
import sys

def main():
    parser = argparse.ArgumentParser(description='Merge protein JSON data while preserving initial structure')
    parser.add_argument('inputs', nargs='+', help='Input protein directories')
    parser.add_argument('--output_dir', help='Directory to save the merged JSON file')
    parser.add_argument('--name', help='Custom name for the output file')
    args = parser.parse_args()

    # Load first input as base structure
    first_input = args.inputs[0]
    base_dir = os.path.dirname(first_input)
    prot = os.path.basename(first_input)
    data_file = os.path.join(base_dir, prot.lower(), f"{prot.lower()}_data.json")

    with open(data_file) as f:
        output = json.load(f)

    # Process subsequent inputs
    for input_path in args.inputs[1:]:
        prot = os.path.basename(input_path)
        base_dir = os.path.dirname(input_path)
        data_file = os.path.join(base_dir, prot.lower(), f"{prot.lower()}_data.json")

        with open(data_file) as f:
            new_data = json.load(f)

        output['sequences'].append(new_data['sequences'][0])

    # Generate sequential protein IDs for all sequences
    for idx, seq in enumerate(output['sequences']):
        seq['protein']['id'] = chr(65 + idx)  # A, B, C...

    # Set metadata
    output['name'] = args.name or '_'.join(os.path.basename(p) for p in args.inputs)
    output['modelSeeds'] = [0]

    # Save result
    os.makedirs(args.output_dir, exist_ok=True)
    output_path = os.path.join(args.output_dir, f"{output['name']}.json")
    
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=4)

if __name__ == '__main__':
    main()