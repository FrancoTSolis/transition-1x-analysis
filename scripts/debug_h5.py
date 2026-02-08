import h5py
import sys

def inspect_h5(path):
    try:
        f = h5py.File(path, 'r')
        print(f"Keys in root: {list(f.keys())}")
        
        if 'data' in f:
            data = f['data']
            print(f"Keys in 'data': {list(data.keys())[:20]}") # Show first 20 formulas
            
            first_formula = list(data.keys())[0]
            print(f"Example reactions for {first_formula}: {list(data[first_formula].keys())[:10]}")
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        inspect_h5(sys.argv[1])
    else:
        print("Usage: python debug_h5.py <h5_file>")
