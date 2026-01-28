# Bioinformatics Assignment: Synthetic Plasmid Assembly Tool

## 1. Project Overview
This project is a bioinformatics tool designed to simulate the automated construction of a plasmid vector. The objective is to take an unknown organism's genome and a specific "design blueprint" to assemble a functional plasmid. 

The tool highlights two key bioinformatics capabilities:
1.  **Algorithmic Discovery:** Identifying a functional Origin of Replication (ORI) from raw genomic data.
2.  **Whitelist Assembly:** Constructing a sequence strictly from requested parts to ensure "negative design" requirements (e.g., removing a specific restriction site) are met.

---

## 2. Methodology & Algorithm

### Part A: Locating the Origin of Replication (ORI)
The tool analyzes the `Input.Fa` (the unknown organism's genome) to find a compatible replication start site using a two-step approach:
* **GC Skew Analysis:** Calculates the cumulative skew ($G-C$) across the genome. In prokaryotes, the global minimum of the GC skew often correlates with the replication origin.
* **Motif Refinement (DnaA Box):** Searches the region around the skew minimum for the consensus DnaA binding sequence (`TTATCCACA`). If found, the ORI is centered on this motif; otherwise, the skew minimum is used.

### Part B: Whitelist Assembly & Verification
The tool parses the `Design.txt` file as a list of instructions. It does not "edit" an existing plasmid; it builds a new one from scratch:
1.  **Database Loading:** Reads sequences for enzymes and markers from `markers.tab`.
2.  **Construction:** Concatenates sequences in the exact order listed in the design file.
3.  **Safety Validation (The "Deletion" Requirement):** A key requirement is determining if a specific site (like EcoRI) was *omitted* from the design. The tool scans the final assembled sequence. If an unrequested `GAATTC` (EcoRI) site appears (e.g., naturally occurring inside a gene or ORI), the tool automatically mutates it (e.g., `GAATTC` $\rightarrow$ `CAATTC`) to ensure the final vector is strictly compliant.

---

## 3. Test Cases

### The "pUC19" Scenario
The assignment provides a specific test case to verify the assembly and deletion logic.

* **Input (`pUC19.fa`):** A raw DNA sequence representing the host organism.
* **Design (`Design_pUC19.txt`):** A blueprint requesting `BamHI`, `HindIII`, `Ampicillin`, and the `ORI`. Crucially, **EcoRI is missing** from this list.
* **Success Criteria:**
    * ✅ **BamHI Check:** The final output must contain `GGATCC` (requested).
    * ✅ **EcoRI Check:** The final output must **NOT** contain `GAATTC` (omitted).

---

## 4. Database Sequences (`markers.tab`)

The tool utilizes a tab-separated database of genetic parts.

### Restriction Enzymes
[cite_start]*Source: Extracted from assignment documentation [cite: 434-481]*

| Enzyme | Sequence | Notes |
| :--- | :--- | :--- |
| **EcoRI** | `GAATTC` | Classic cloning site (Deleted in Test Case) |
| **BamHI** | `GGATCC` | 5' overhangs |
| **HindIII** | `AAGCTT` | Standard MCS site |
| **PstI** | `CTGCAG` | 3' overhangs |
| **KpnI** | `GGTACC` | 3' overhangs |
| **SacI** | `GAGCTC` | pUC-type MCS |
| **SalI** | `GTCGAC` | Promoter/ORF junctions |
| **XbaI** | `TCTAGA` | Modular cloning |
| **NotI** | `GCGGCCGC` | Rare 8-bp cutter |
| **SmaI** | `CCCGGG` | Blunt ends |
| **BsaI** | `GGTCTC` | Type IIS (Golden Gate) |
| **BbsI** | `GAAGAC` | Type IIS |
| **BsmBI** | `CGTCTC` | Type IIS |

### Selectable Markers & Genes
[cite_start]*Source: Standard biological reference sequences corresponding to assignment descriptions [cite: 482-501]*

| Marker | Role |
| :--- | :--- |
| **AmpR (bla)** | Ampicillin resistance ($\beta$-lactamase). Standard pUC marker. |
| **KanR (nptII)** | Kanamycin resistance. Robust selection. |
| **CmR (cat)** | Chloramphenicol resistance. Low-copy plasmids. |
| **TetR (tetA)** | Tetracycline resistance. pBR322 lineage. |
| **SpecR (aadA)** | Spectinomycin resistance. Broad-host-range. |
| **lacZα** | Blue/White screening fragment. |
