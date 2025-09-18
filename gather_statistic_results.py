import os
import glob

EXPORT_DIR = "./Exports"
COLLECT_DIR = "./CollectedAnalytics"
ARCHS = ["cuda-8", "mpi-8", "single-core"]

os.makedirs(COLLECT_DIR, exist_ok=True)

for arch in ARCHS:
    arch_path = os.path.join(COLLECT_DIR, arch)
    os.makedirs(arch_path, exist_ok=True)

    case_paths = glob.glob(os.path.join(EXPORT_DIR, arch, "*"))

    for case_path in case_paths:
        case_name = os.path.basename(case_path)
        stats_path = os.path.join(case_path, "SimulationStats")
        if not os.path.isdir(stats_path):
            continue

        rank_files = glob.glob(os.path.join(stats_path, "statistics_rank*.dat"))
        if not rank_files:
            continue

        output_file = os.path.join(arch_path, f"{case_name}.dat")

        # Initialize variables
        init_times = []
        run_times = []
        element_sizes = []
        final_time = None
        timestep = None
        num_elements = 0
        num_local_dofs = 0
        operator_total_local = 0
        operator_total_local_ghost = 0
        operator_nonzero = 0

        # Single-core: just read rank0
        if arch == "single-core":
            rank_file = rank_files[0]
            with open(rank_file, "r") as f:
                for line in f:
                    if ":" not in line:
                        continue
                    key, value = line.split(":", 1)
                    value = value.strip().split()[0]
                    try:
                        if key.startswith("Initialization Time"):
                            init_times.append(float(value))
                        elif key.startswith("Simulation Run Time"):
                            run_times.append(float(value))
                        elif key.startswith("Final Time"):
                            final_time = float(value)
                        elif key.startswith("Time Step"):
                            timestep = float(value)
                        elif key.startswith("Number of Mesh Elements"):
                            num_elements += int(value)
                        elif key.startswith("Average Element Size in Mesh"):
                            element_sizes.append(float(value))
                        elif key.startswith("Number of Local Degrees of Freedom"):
                            num_local_dofs += int(value)
                        elif key.startswith("Operator Total Elements for Local Degrees of Freedom"):
                            operator_total_local = int(value)
                        elif key.startswith("Operator Total Elements for Local and Ghost Degrees of Freedom"):
                            operator_total_local_ghost = int(value)
                        elif key.startswith("Number of Operator Non-Zero Elements"):
                            operator_nonzero = int(value)
                    except ValueError:
                        pass

        # Multi-rank: read all rank files and aggregate
        else:
            for rank_file in rank_files:
                with open(rank_file, "r") as f:
                    for line in f:
                        if ":" not in line:
                            continue
                        key, value = line.split(":", 1)
                        value = value.strip().split()[0]
                        try:
                            if key.startswith("Initialization Time"):
                                init_times.append(float(value))
                            elif key.startswith("Simulation Run Time"):
                                run_times.append(float(value))
                            elif key.startswith("Final Time") and final_time is None:
                                final_time = float(value)  # same for all ranks
                            elif key.startswith("Time Step") and timestep is None:
                                timestep = float(value)  # same for all ranks
                            elif key.startswith("Number of Mesh Elements"):
                                num_elements += int(value)
                            elif key.startswith("Average Element Size in Mesh"):
                                element_sizes.append(float(value))
                            elif key.startswith("Number of Local Degrees of Freedom"):
                                num_local_dofs += int(value)
                            elif key.startswith("Operator Total Elements for Local Degrees of Freedom"):
                                operator_total_local += int(value)
                            elif key.startswith("Operator Total Elements for Local and Ghost Degrees of Freedom"):
                                operator_total_local_ghost += int(value)
                            elif key.startswith("Number of Operator Non-Zero Elements"):
                                operator_nonzero += int(value)
                        except ValueError:
                            pass

        # Compute averages where needed
        avg_init_time = sum(init_times) / len(init_times) if init_times else 0.0
        avg_run_time = sum(run_times) / len(run_times) if run_times else 0.0
        avg_element_size = sum(element_sizes) / len(element_sizes) if element_sizes else 0.0

        # Write aggregated data
        with open(output_file, "w") as out:
            out.write(f"Initialization Time: {avg_init_time:.6e} (s)\n")
            out.write(f"Simulation Run Time: {avg_run_time:.6e} (s)\n")
            out.write(f"Final Time: {final_time} (ns)\n")
            out.write(f"Time Step: {timestep} (ns)\n")
            out.write(f"Number of Mesh Elements: {num_elements}\n")
            out.write(f"Average Element Size in Mesh: {avg_element_size}\n")
            out.write(f"Number of Local Degrees of Freedom: {num_local_dofs}\n")
            out.write(f"Operator Total Elements for Local Degrees of Freedom: {operator_total_local}\n")
            out.write(f"Operator Total Elements for Local and Ghost Degrees of Freedom: {operator_total_local_ghost}\n")
            out.write(f"Number of Operator Non-Zero Elements: {operator_nonzero}\n")
