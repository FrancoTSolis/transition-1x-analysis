# QuantumATK Quick Start

## Installation Complete! ✓

QuantumATK X-2025.06 is installed at: `/media/francosolis/newdrive/qatk/`

## Initial Setup (One-time only)

1. Source the updated shell configuration:
   ```bash
   source ~/.zshrc
   ```

2. Start the persistent license tunnel:
   ```bash
   qatk-tunnel-start
   ```

That's it! The tunnel will persist across shell sessions.

## Daily Usage

After the initial setup, you can use QuantumATK normally in any shell:

```bash
# Interactive Python
atkpython
# or simply
qatk

# Run a script
atkpython my_script.py

# Use GPU version
qatk-gpu my_script.py

# Use MPI version
qatk-mpi my_script.py

# Start GUI
qatk-gui
```

## Tunnel Management

```bash
# Check tunnel status
qatk-tunnel-status

# Stop tunnel (rarely needed)
qatk-tunnel-stop

# Restart tunnel
qatk-tunnel-start
```

## Available Aliases

- `qatk` → `atkpython`
- `qatk-gpu` → `atkpython_gpu`
- `qatk-mpi` → `atkpython_system-mpi`
- `qatk-gui` → `quantumatk`

## Documentation

See `QuantumATK_Setup_Manual.md` for detailed documentation.

## License Server

- Server: `pyrito.cs.ucla.edu:27020`
- Access: Through persistent SSH tunnel via `scai1.cs.ucla.edu`
- Local configuration: `SNPSLMD_LICENSE_FILE=27020@localhost`
