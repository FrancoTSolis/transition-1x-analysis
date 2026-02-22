**Checklist for Collaboration with Prof. Martin Head-Gordon's Group** 

**Subject:** Request for CCSD(T) and MP2 Implementation Examples for Transition State Benchmarks 

**Context:** We are benchmarking quantum circuit initialization strategies (LUCJ ansatz) and need high-precision classical amplitudes. We require "jump-start" support to run efficient single-point energy calculations on transition state geometries derived from literature. 

**Primary Request:** We need working input examples and software guidance for the following two specific regimes: 

**1. Gold Standard: CCSD(T) for Medium Molecules (~20 Atoms)** 

We need a **running example (input deck & submission script)** for a ~20-atom transition state molecule using **CCSD(T)** with the following specifications: 

- **Basis Set 1:** 6-31G(d) 

- **Basis Set 2:** def2-TZVPD 

- **Requirements:** 

- Optimal memory/core settings for standard HPC nodes (e.g., 256GB RAM). 

- Recommended software flags to ensure convergence on transition states. 

- Estimated wall-time for these specific basis sets. 

**2. Silver Tier: MP2 Approximation for Larger Molecules (~30 Atoms)** 

We need a **running example (input deck & submission script)** for a ~30-atom transition state molecule using **MP2** as a cost-effective proxy for CCSD(T). 

- **Goal:** Use MP2 amplitudes to approximate the coupled-cluster initialization. 

- **Requirements:** 

- Input setup that aligns with the CCSD(T) reference (same frozen core settings, etc.). 

- Guidance on validity: Confirmation of MP2's reliability for these transition metal/organic transition states compared to CCSD(T). 

**3. Software Environment & Execution** 

- **Preferred Code:** Please confirm if we should use Q-Chem (or another package preferred by your group). 

- **Environment:** If using Q-Chem, do you have a recommended module/container for SDSC Expanse? 

- **Output Handling:** We specifically need to extract **t_1 and t_2 amplitudes** and **molecular orbital coefficients**. Please point to the relevant output sections or keywords to print these explicitly. 

 

 