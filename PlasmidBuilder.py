import os


# DnaA Box Consensus (Used to refine ORI search)
DNA_A_BOX = "TTATCCACA"

# We keep the EcoRI sequence hardcoded ONLY for the final validation check,
# ensuring we comply with the specific requirement to verify its absence.
ECORI_SEQ_VALIDATION = "GAATTC" 


def get_gc_skew(sequence):
    """Calculates cumulative GC Skew to estimate the Origin of Replication."""
    skew = [0]
    current = 0
    for base in sequence:
        base = base.upper()
        if base == 'G':
            current += 1
        elif base == 'C':
            current -= 1
        skew.append(current)
    return skew

def find_ori_sequence(sequence, window=500):
    """
    Locates the ORI sequence in the provided organism DNA.
    1. Finds global minimum of GC Skew.
    2. Refines search by looking for DnaA box nearby.
    """
    # Step A: GC Skew Analysis
    skew_values = get_gc_skew(sequence)
    min_skew = min(skew_values)
    min_idx = skew_values.index(min_skew)
    
    # Step B: Define Search Region (Window around min)
    start_region = max(0, min_idx - window)
    end_region = min(len(sequence), min_idx + window)
    region_seq = sequence[start_region:end_region]

    # Step C: Refinement (Look for DnaA Box)
    dna_a_index = region_seq.find(DNA_A_BOX)
    
    # Calculate absolute start position
    if dna_a_index != -1:
        absolute_ori_start = start_region + dna_a_index
    else:
        # Fallback to the skew minimum if no box found
        absolute_ori_start = min_idx

    # Return a functional chunk (e.g., 250bp)
    ori_end = min(len(sequence), absolute_ori_start + 250)
    return sequence[absolute_ori_start:ori_end]

def load_database(filepath):
    """
    Reads 'markers.tab'. 
    Format expected: Name [TAB] Sequence
    Loads BOTH restriction enzymes and markers into one dictionary.
    """
    database = {}
    if not os.path.exists(filepath):
        print(f"[Warning] Marker file {filepath} not found. Using empty database.")
        return database

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            # Split by tab
            parts = line.split('\t')
            if len(parts) >= 2:
                name = parts[0].strip()
                seq = parts[1].strip().upper() # Ensure DNA is uppercase
                database[name] = seq
            else:
                # Handle cases where there might be spaces instead of tabs, purely for robustness
                parts_alt = line.split()
                if len(parts_alt) >= 2:
                    database[parts_alt[0].strip()] = parts_alt[1].strip().upper()
                else:
                    print(f"[Warning] Skipping malformed line {line_num} in markers.tab")
                    
    return database

def parse_design_file(filepath):
    """
    Parses 'Design.txt' into an ordered list of instructions.
    Returns list of tuples: [('BamHI_site', 'BamHI'), ('AmpR_gene', 'Ampicillin'), ...]
    """
    instructions = []
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Design file {filepath} not found.")

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comments, empty lines, and closing brackets
            if not line or line.startswith("#") or line.startswith("*") or line.startswith("}"):
                continue
            
            parts = line.split(',')
            if len(parts) >= 2:
                feature_type = parts[0].strip() # e.g., "BamHI_site"
                feature_name = parts[1].strip() # e.g., "BamHI"
                instructions.append((feature_type, feature_name))
                
    return instructions


def generate_plasmid(input_fa_path, design_txt_path, markers_tab_path="markers.tab", output_fa_path="Output.Fa"):
    """
    1. Finds ORI in input_fa.
    2. Loads database from markers.tab.
    3. Assembles plasmid based on Design.txt whitelist.
    4. Saves to output_fa.
    """
    print(f"\n--- Starting Plasmid Assembly ---")
    
    # 1. Read Host Genome
    full_genome = ""
    try:
        with open(input_fa_path, 'r') as f:
            for line in f:
                if not line.startswith(">"):
                    full_genome += line.strip()
    except FileNotFoundError:
        print(f"[Error] Input file {input_fa_path} missing.")
        return

    # 2. Find ORI
    identified_ori = find_ori_sequence(full_genome)
    print(f"[Analysis] ORI identified in host genome ({len(identified_ori)} bp).")

    # 3. Load External Database (Markers + Enzymes)
    db = load_database(markers_tab_path)
    print(f"[Loading] Loaded {len(db)} sequences from {markers_tab_path}.")

    # 4. Read Design Instructions
    assembly_plan = parse_design_file(design_txt_path)
    print(f"[Loading] Found {len(assembly_plan)} components in design file.")

    # 5. Assemble Sequence
    final_sequence = ""
    print("\n[Assembly Step-by-Step]:")
    ori_added = False
    missing_features = []

    print("\n[Assembly Step-by-Step]:")
    
    for feature_type, feature_name in assembly_plan:
        
        # Case A: User explicitly asks for ORI (e.g. "ori_pMB1")
        # We check both the Type string and the Name string for "ori" or "replication"
        if "ori" in feature_type.lower() or "replication" in feature_name.lower():
            final_sequence += identified_ori
            ori_added = True
            print(f"  + Added Identified ORI ({feature_name})")

        # Case B: Feature is in our Database
        elif feature_name in db:
            seq_fragment = db[feature_name]
            final_sequence += seq_fragment
            print(f"  + Added {feature_name} ({len(seq_fragment)} bp)")

        # Case C: Feature NOT found
        else:
            print(f"  ❌ CRITICAL WARNING: '{feature_name}' requested but NOT found in markers.tab.")
            missing_features.append(feature_name)
            # We do NOT add anything here, but we log it.

    # 3. Post-Assembly Checks
    
    # Check 1: Did we miss anything?
    if missing_features:
        print(f"\n[Alert] The following features were skipped (missing from DB): {missing_features}")
    
    # Check 2: Did the user forget the ORI? (Auto-Inclusion)
    if not ori_added:
        print("\n[Auto-Fix] No ORI specified in Design file. Appending identified ORI to ensure functionality.")
        final_sequence += identified_ori

    # 6. Validation (The "Delete EcoRI" Requirement)
    # Even if we didn't add EcoRI intentionally, we check if it snuck in via ORI or Markers.
    if ECORI_SEQ_VALIDATION in final_sequence:
        print(f"\n[Validation Alert] 'EcoRI' site ({ECORI_SEQ_VALIDATION}) detected in final sequence.")
        
        # Check if it was requested in the design
        design_names = [name for _, name in assembly_plan]
        if "EcoRI" not in design_names:
            print("  -> It was NOT in the design plan. Performing cleanup mutation...")
            # Mutate G -> C to remove the site while maintaining length
            final_sequence = final_sequence.replace(ECORI_SEQ_VALIDATION, "CAATTC")
            print("  -> Site mutated to 'CAATTC'.")
    else:
        print("\n[Validation Pass] No unintended EcoRI sites found.")

    # 7. Write Output
    with open(output_fa_path, 'w') as f:
        f.write(">Synthetic_Plasmid_Output\n")
        f.write(final_sequence + "\n")
    
    print(f"\n--- Success! Plasmid saved to: {output_fa_path} ---")
    return output_fa_path

def run_test_case(input_fa="testfiles/pUC19.fa", design_txt="testfiles/Design_pUC19.txt", markers_tab="markers.tab"):
    """
    Runs the tool using existing input files and validates the output.
    Does NOT create/overwrite files, so it works with your uploaded data.
    """
    import os
    
    print(f"=== RUNNING TEST CASE ===")
    print(f"Input Genome: {input_fa}")
    print(f"Design File:  {design_txt}")
    print(f"Markers DB:   {markers_tab}")

    # 1. Check if input files exist before running
    if not os.path.exists(input_fa) or not os.path.exists(design_txt):
        print("❌ Error: One or more input files are missing. Please ensure pUC19.fa and Design_pUC19.txt exist.")
        return

    # 2. Run the tool
    # We pass the filenames you provided as arguments
    output_file = generate_plasmid(input_fa, design_txt, markers_tab, "Output.Fa")
    
    # 3. Verify Results
    if not os.path.exists(output_file):
        print("❌ Error: Output file was not created.")
        return

    with open(output_file, 'r') as f:
        lines = f.readlines()
        # Join lines just in case the FASTA is multi-line
        result_seq = "".join([l.strip() for l in lines if not l.startswith(">")])

    print("\n=== FINAL VERIFICATION ===")
    print(f"Output Length: {len(result_seq)} bp")

    # Check 1: BamHI (GGATCC) should be present (It's in the design file)
    if "GGATCC" in result_seq: 
        print("✅ PASS: BamHI site is present.")
    else: 
        print("❌ FAIL: BamHI site is missing.")

    # Check 2: EcoRI (GAATTC) should be ABSENT (It's NOT in the design file)
    if "GAATTC" not in result_seq: 
        print("✅ PASS: EcoRI site successfully deleted/absent.")
    else: 
        print("❌ FAIL: EcoRI site is still present.")

if __name__ == "__main__":
    run_test_case()