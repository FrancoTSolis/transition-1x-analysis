# QuantumATK Installation and Usage Manual

## Overview

This manual describes how to use QuantumATK on a local machine when the license server (`pyrito.cs.ucla.edu`) is only accessible from a whitelisted server (`scai1.cs.ucla.edu`). The solution uses SSH tunneling to forward license requests through the whitelisted server.

## Prerequisites

1. **QuantumATK Installation**: Installed at `/media/francosolis/newdrive/qatk/`
2. **SSH Access**: Valid credentials for `fts@scai1.cs.ucla.edu`
3. **License Server**: `27020@pyrito.cs.ucla.edu` (accessible only from scai1)

## Quick Start

### Initial Setup (One-time only)

1. **Update your shell configuration**:
   ```bash
   source ~/.zshrc
   ```

2. **Start the persistent license tunnel**:
   ```bash
   qatk-tunnel-start
   ```

### Daily Usage

After the initial setup, you can use QuantumATK directly in any shell:
```bash
atkpython    # or use the alias: qatk
```

The tunnel runs persistently in the background and survives shell sessions.

## Available QuantumATK Binaries

The following binaries are available in `/media/francosolis/newdrive/qatk/tools/quantumatk/X-2025.06/bin/`:

- **atkpython**: Main QuantumATK Python interpreter
- **atkpython_gpu**: GPU-accelerated version
- **atkpython_system-mpi**: MPI-enabled version for parallel computing
- **quantumatk**: QuantumATK GUI application

## Usage Examples

### Basic Usage
```bash
# Start interactive atkpython
atkpython
# or use the alias
qatk

# Run a Python script
atkpython script.py

# Run with specific Python commands
atkpython -c "print('Hello from QuantumATK')"
```

### Using Different Binaries
```bash
# Use GPU version
atkpython_gpu script.py
# or use the alias
qatk-gpu script.py

# Use MPI version
atkpython_system-mpi script.py
# or use the alias
qatk-mpi script.py

# Start the GUI
quantumatk
# or use the alias
qatk-gui
```

### Environment Variables

The following environment variables are automatically set in your `.zshrc`:

```bash
# QuantumATK installation directory
QATK_HOME="/media/francosolis/newdrive/qatk/tools/quantumatk/X-2025.06"

# QuantumATK binaries in PATH
PATH="$QATK_HOME/bin:$PATH"

# License server configuration (using localhost tunnel)
SNPSLMD_LICENSE_FILE="27020@localhost"
```

## How It Works

### The License Server Problem
- The Synopsys license server at `pyrito.cs.ucla.edu` only accepts connections from whitelisted IP addresses
- Your local machine is not on the whitelist, but `scai1.cs.ucla.edu` is
- Direct connection attempts from your local machine will be rejected

### The SSH Tunnel Solution
1. SSH creates encrypted tunnels from your local machine to `scai1.cs.ucla.edu`
2. Local ports (27020-27025, 49200-49204) are forwarded through these tunnels
3. When QuantumATK connects to `localhost:27020`, the connection is forwarded through scai1
4. The license server sees the connection coming from scai1 (whitelisted) and accepts it

### Why Not SOCKS Proxy?
- SOCKS proxy works at a different network layer and is great for HTTP/HTTPS traffic
- FlexNet licensing (used by Synopsys) requires direct TCP connections on specific ports
- FlexNet may also use UDP and dynamic port allocation, which SOCKS doesn't handle well
- SSH port forwarding creates proper TCP tunnels that FlexNet can use transparently

## Troubleshooting

### SSH Tunnel Issues
```bash
# Check if tunnel is running
ps aux | grep -E "ssh.*27020.*pyrito"

# Kill existing tunnels and restart
pkill -f "ssh.*27020.*pyrito"
./setup_ssh_tunnel_flexnet.sh
```

### License Connection Test
```bash
# Test connection to license server
nc -zv localhost 27020
nc -zv localhost 27021

# These should show "Connection succeeded"
```

### Environment Variables
```bash
# Check current license configuration
echo $SNPSLMD_LICENSE_FILE

# Should output: 27020@localhost
```

### Common Issues

1. **"Cannot connect to license server"**
   - Ensure SSH tunnel is running
   - Check that SNPSLMD_LICENSE_FILE is set to `27020@localhost`
   - Verify SSH access to scai1 is working

2. **"Operation now in progress" error**
   - This typically means the connection is timing out
   - Make sure you're using SSH tunneling, not SOCKS proxy

3. **Multiple license server entries**
   - Check for conflicting environment variables
   - Look for `.quantumatk/licenses/snps.licconf` in your home directory
   - Ensure only one license server is configured

## Files and Scripts

### Essential Files
- `qatk_tunnel_daemon.sh`: Manages the persistent SSH tunnel
  - `qatk-tunnel-start`: Start the tunnel
  - `qatk-tunnel-stop`: Stop the tunnel
  - `qatk-tunnel-status`: Check tunnel status

### Configuration Locations
- QuantumATK installation: `/media/francosolis/newdrive/qatk/tools/quantumatk/X-2025.06/`
- User configuration: `~/.quantumatk/licenses/snps.licconf`
- Shell configuration: `~/.zshrc` (contains PATH and aliases)
- Tunnel PID file: `~/.quantumatk/tunnel.pid`
- Tunnel log file: `~/.quantumatk/tunnel.log`

## Best Practices

1. **One-time setup**: Run `qatk-tunnel-start` once after boot/login
2. **Keep tunnel alive**: The tunnel daemon runs persistently in the background
3. **One tunnel only**: The daemon script prevents multiple tunnels
4. **Use aliases**: The shell aliases make commands shorter and easier to remember

## Advanced Usage

### Running on Multiple Nodes
For MPI jobs across multiple nodes, each node needs access to the license server:
```bash
# On each compute node
./setup_ssh_tunnel_flexnet.sh
export SNPSLMD_LICENSE_FILE=27020@localhost

# Then run your MPI job
mpirun -np 4 atkpython_system-mpi script.py
```

### Persistent Tunnel with systemd
For a persistent tunnel that survives reboots, create a systemd service (requires root):
```bash
sudo tee /etc/systemd/system/qatk-license-tunnel.service << EOF
[Unit]
Description=QuantumATK License Tunnel
After=network.target

[Service]
Type=simple
User=$USER
ExecStart=/usr/bin/ssh -N -L 27020:pyrito.cs.ucla.edu:27020 -L 27021:pyrito.cs.ucla.edu:27021 fts@scai1.cs.ucla.edu
Restart=always
RestartSec=30

[Install]
WantedBy=multi-user.target
EOF

sudo systemctl enable qatk-license-tunnel.service
sudo systemctl start qatk-license-tunnel.service
```

## Contact and Support

For issues specific to:
- **SSH/Network access**: Contact your system administrator
- **License server**: Contact your license administrator
- **QuantumATK software**: Visit https://forum.quantumatk.com or https://solvnet.synopsys.com

## Sample .zshrc setup

```
# QuantumATK Configuration
# ========================

# QuantumATK Installation Path
export QATK_HOME="/media/francosolis/newdrive/qatk/tools/quantumatk/X-2025.06"

# Add QuantumATK binaries to PATH
export PATH="$QATK_HOME/bin:$PATH"

# License configuration (using localhost when tunnel is active)
export SNPSLMD_LICENSE_FILE="27020@localhost"

# Function to ensure QuantumATK tunnel is running
qatk_ensure_tunnel() {
    local tunnel_script="/media/francosolis/newdrive/to_fang_20250826/qatk_tunnel_daemon.sh"

    if [ -f "$tunnel_script" ]; then
        if ! "$tunnel_script" status >/dev/null 2>&1; then
            echo "Starting QuantumATK license tunnel..."
            "$tunnel_script" start
        fi
    else
        echo "Warning: QuantumATK tunnel script not found at $tunnel_script"
    fi
}

# Aliases for QuantumATK commands
alias qatk='atkpython'
alias qatk-gpu='atkpython_gpu'
alias qatk-mpi='atkpython_system-mpi'
alias qatk-gui='quantumatk'
alias qatk-tunnel-start='/media/francosolis/newdrive/to_fang_20250826/qatk_tunnel_daemon.sh start'
alias qatk-tunnel-stop='/media/francosolis/newdrive/to_fang_20250826/qatk_tunnel_daemon.sh stop'
alias qatk-tunnel-status='/media/francosolis/newdrive/to_fang_20250826/qatk_tunnel_daemon.sh status'
```

---
Last updated: September 30, 2025
