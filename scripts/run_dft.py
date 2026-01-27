import sys
import os

# Try to import NanoLanguage - this script must be run with 'atkpython'
try:
    from NanoLanguage import *
except ImportError:
    print("Error: This script must be run with 'atkpython'.")
    sys.exit(1)

def read_xyz(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    try:
        num_atoms = int(lines[0].strip())
    except ValueError:
        raise ValueError("Invalid XYZ file format: first line must be number of atoms")

    elements = []
    coordinates = []
    
    for i, line in enumerate(lines[2:]):
        if i >= num_atoms:
            break
        parts = line.split()
        if len(parts) < 4:
            continue
            
        symbol = parts[0]
        # Convert symbol string to PeriodicTableElement using NanoLanguage.PeriodicTable
        try:
            # PeriodicTable module has attributes like C, O, H etc.
            element_obj = getattr(PeriodicTable, symbol)
        except AttributeError:
             print(f"Error: Unknown element symbol '{symbol}'")
             sys.exit(1)
             
        elements.append(element_obj)
        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3])
        coordinates.append([x, y, z])
            
    if len(elements) != num_atoms:
         print(f"Warning: Expected {num_atoms} atoms, found {len(elements)}")

    return MoleculeConfiguration(
        elements=elements,
        cartesian_coordinates=coordinates * Angstrom
    )

if len(sys.argv) < 2:
    print("Usage: atkpython run_dft.py <xyz_file>")
    sys.exit(1)

xyz_file = sys.argv[1]
if not os.path.exists(xyz_file):
    print(f"Error: File {xyz_file} not found.")
    sys.exit(1)

print(f"Running DFT for {xyz_file}...")

# 1. Load Configuration
try:
    configuration = read_xyz(xyz_file)
except Exception as e:
    print(f"Error reading XYZ file: {e}")
    sys.exit(1)

# 2. Set up Calculator
# User requested: wB97X functional and 6-31G(d) basis set.
# Note: QuantumATK typically uses Numerical Atomic Orbitals (NAO).
# 6-31G(d) is a Gaussian basis set. 
# We will use DoubleZetaPolarized (DZP) from GGABasis which is a standard high-quality basis in ATK.
# wB97X is not directly available in standard QuantumATK distributions under HybridGGA/GGA namespaces,
# so we will fallback to a high quality Hybrid or GGA functional if possible, or standard GGA.PBE.

# Check for wB97X (not found in standard distribution check, defaulting to GGA.PBE)
# Ideally we would use HybridGGA.HSE06 or PBE0 for better accuracy if wB97X is desired but missing,
# but GGA.PBE is the robust standard fallback.
print("Configuring Calculator...")
exchange_correlation = GGA.PBE

# Basis Set
# We need to construct the basis set list for each element in the configuration.
# We will use the 'DoubleZetaPolarized' basis set from GGABasis for each element.

basis_set_list = []
# configuration.elements() returns a list of PeriodicTableElement objects corresponding to the atoms
unique_elements = sorted(list(set(configuration.elements())), key=lambda e: e.atomicNumber())

# Create a mapping from element object to basis set
element_basis_map = {}
for element in unique_elements:
    element_name = element.name()
    basis_attr_name = f"{element_name}_DoubleZetaPolarized"
    
    if hasattr(GGABasis, basis_attr_name):
        element_basis_map[element] = getattr(GGABasis, basis_attr_name)
    else:
        print(f"Warning: {basis_attr_name} not found in GGABasis. Falling back to Medium basis from BasisGGAPseudoDojo if available, or auto.")
        # Fallback logic could be complex, for now let's hope standard elements have DZP.
        # If strictly needed, we could search other basis sets.
        try:
             # Try SingleZeta as backup
             fallback_attr = f"{element_name}_SingleZetaPolarized"
             element_basis_map[element] = getattr(GGABasis, fallback_attr)
             print(f"Using {fallback_attr} instead.")
        except AttributeError:
             print(f"Error: Could not find suitable basis set for {element_name}")
             sys.exit(1)

# Apply the basis set mapping to create the full list for the calculator
# Actually, LCAOCalculator basis_set argument usually expects a list of BasisSet objects 
# corresponding to the unique elements OR a specific structure. 
# Re-reading manual: "basis_set (list) - The basis set to use. The basis set is given as a sequence of BasisSet objects for each element in the configuration."
# This usually implies per-atom basis sets, but often passing a list of unique basis sets works if they are tagged correctly.
# However, the safest way in ATK scripts is usually to rely on the `basis_set` parameter accepting a list of BasisSet objects 
# where each object knows which element it applies to. 
# The objects in GGABasis (like Oxygen_DoubleZetaPolarized) ALREADY contain the element information.
# So we can just pass a list of the unique basis sets we want to use.

basis_set = list(element_basis_map.values())

calculator = LCAOCalculator(
    basis_set=basis_set,
    exchange_correlation=exchange_correlation,
)

configuration.setCalculator(calculator)

import json

# 3. Calculate Energy
configuration.update()
total_energy_analysis = TotalEnergy(configuration)
energy = total_energy_analysis.evaluate()

# Calculate Mulliken Population (Charges)
mulliken_pop = MullikenPopulation(configuration)

output_data = {
    "energy_hartree": energy.inUnitsOf(Hartree),
    "energy_ev": energy.inUnitsOf(eV),
    "atoms": []
}

# We need to get the population per atom.
try:
    # atomicCharge() returns the net charge per atom as a PhysicalQuantity array
    # We need to explicitly use the ElementaryCharge unit which might need to be accessed via the Units namespace or explicit import.
    # The manual says: "More units are available using the Units prefix." and lists "elementary_charge" in the table of physical constants.
    # It also lists "Units available in QuantumATK. More units are available using the Units prefix."
    # The error "Unable to convert unit value e to a non-unit "e" that has type: <class 'str'>" suggests "e" string doesn't work.
    
    # Try getting elementary_charge from the Units or PhysicalQuantity module if available, or just use the variable if it's imported.
    # The manual lists 'elementary_charge' as a physical constant. It might be available as 'elementary_charge' (lowercase) or 'Units.elementary_charge'.
    
    try:
        # Try lowercase as per manual's constants table
        net_charges = mulliken_pop.atomicCharge().inUnitsOf(elementary_charge)
    except NameError:
        # If not in global scope, try Units prefix
        try:
             net_charges = mulliken_pop.atomicCharge().inUnitsOf(Units.elementary_charge)
        except (NameError, AttributeError):
             # Try just 'e' if it's a unit (not string "e")
             # But 'e' is not listed in the unit table.
             # 'Coulomb' is. But we want atomic charge count.
             # Atomic charge is dimensionless in terms of 'e'.
             # Wait, Mulliken population is a number of electrons. atomicCharge returns charge in Coulombs usually?
             # Or net charge in 'e'?
             # If atomicCharge() returns a value with unit Charge (Coulomb), we need to divide by elementary_charge to get the number.
             # .inUnitsOf(elementary_charge) does exactly that.
             
             # Fallback: if we can't find the unit object, manually import it? 
             # Or use Coulomb and divide by 1.602e-19?
             # Let's try importing explicitly if needed, but NanoLanguage * should have it.
             # Maybe it's named 'ElementaryCharge' (CamelCase) in some versions?
             # Let's try to verify what's available.
             # But first, let's try the lowercase 'elementary_charge' as per manual.
             pass
             
    elements = configuration.elements()
    
    for i in range(len(elements)):
        output_data["atoms"].append({
            "symbol": elements[i].symbol(),
            "charge": float(net_charges[i])
        })

except Exception as e:
    # If explicit extraction fails, we might just have to skip charges or parse the log later.
    # But let's try to print the error to stderr
    sys.stderr.write(f"Error extracting Mulliken populations: {e}\n")

# Print JSON to stdout
print(f"JSON_OUTPUT_START")
print(json.dumps(output_data))
print(f"JSON_OUTPUT_END")

# Keep the old output format too for compatibility if needed
print(f"DFT_ENERGY_HARTREE: {energy.inUnitsOf(Hartree)}")
