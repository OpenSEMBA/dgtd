import os
import shutil

# Root paths
exports_root = "Exports/mpi-8"
collected_root = "CollectedAnalytics/Probes"

# Ensure collected folder exists
os.makedirs(collected_root, exist_ok=True)

# Header to prepend to PointProbe0.dat
header = """PointProbe ID 0
Spatial Position (X, Y, Z) 
0.250000 0.750000 0.000000
Time (s) // Ex // Ey // Ez // Hx // Hy // Hz 
"""

# Sweep through all casename folders
for casename in os.listdir(exports_root):
    case_path = os.path.join(exports_root, casename, "PointProbes")
    if not os.path.isdir(case_path):
        continue  # skip non-folders or missing PointProbes

    for probe_id in [0, 1]:
        probe_file = f"PointProbe{probe_id}.dat"
        src_path = os.path.join(case_path, probe_file)

        if not os.path.exists(src_path):
            continue  # skip if probe missing

        # Fix PointProbe0.dat
        if probe_id == 0:
            with open(src_path, "r") as f:
                content = f.read()
            # Check if header is already added
            if not content.startswith("PointProbe ID 0"):
                content = header + content
                with open(src_path, "w") as f:
                    f.write(content)

        # Destination file with casename prefix
        dst_name = f"{casename}_{probe_file}"
        dst_path = os.path.join(collected_root, dst_name)

        shutil.copy2(src_path, dst_path)
        print(f"Collected {dst_name}")
