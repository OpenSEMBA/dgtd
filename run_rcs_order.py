import json
import os
import subprocess

# Base JSON file
BASE_JSON = "3D_RCS_Sphere.json"

# Output directory
OUTPUT_DIR = "./testData/maxwellInputs"

# Program path
PROGRAM = "./build/gnu-release/bin/opensemba_dgtd"

def modify_json(order_value):
    order_suffix = f"_O{order_value}"
    new_filename = f"3D_RCS_Sphere{order_suffix}.json"
    new_folder = os.path.join(OUTPUT_DIR, f"3D_RCS_Sphere{order_suffix}")
    
    # Ensure the output directory exists
    os.makedirs(new_folder, exist_ok=True)
    
    # Read base JSON
    with open(BASE_JSON, 'r') as f:
        data = json.load(f)
    
    # Modify JSON fields
    data['solver_options']['order'] = order_value
    data['model']['filename'] = data['model']['filename'].replace("3D_RCS_Sphere", f"3D_RCS_Sphere{order_suffix}")
    data['probes']['farfield'][0]['name'] = data['probes']['farfield'][0]['name'].replace("circle_1m", f"circle_1m{order_suffix}")
    
    # Save modified JSON
    new_json_path = os.path.join(new_folder, new_filename)
    with open(new_json_path, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"Created: {new_json_path}")
    return new_json_path

def run_program(json_path):
    if os.path.isfile(json_path):
        print(f"Running program with input: {json_path}")
        process = subprocess.run([PROGRAM, "-i", json_path])
        
        if process.returncode != 0:
            print(f"Error: Program failed with exit code {process.returncode}")
            exit(1)
    else:
        print(f"Error: JSON file not found at {json_path}")
        exit(1)

# Iterate over order values 1, 2, 3
for order_value in [1, 2, 3]:
    json_path = modify_json(order_value)
    run_program(json_path)
