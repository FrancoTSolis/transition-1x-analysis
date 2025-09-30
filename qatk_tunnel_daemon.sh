#!/bin/bash

# QuantumATK License Tunnel Daemon
# This script creates a persistent SSH tunnel that survives shell sessions

TUNNEL_PID_FILE="$HOME/.quantumatk/tunnel.pid"
TUNNEL_LOG_FILE="$HOME/.quantumatk/tunnel.log"

# Create directory if it doesn't exist
mkdir -p "$HOME/.quantumatk"

start_tunnel() {
    # Check if tunnel is already running
    if is_tunnel_running; then
        echo "Tunnel is already running (PID: $(cat "$TUNNEL_PID_FILE"))"
        return 0
    fi
    
    echo "Starting QuantumATK license tunnel..."
    
    # Start SSH tunnel in background with automatic restart
    nohup ssh -N \
        -L 27020:pyrito.cs.ucla.edu:27020 \
        -L 27021:pyrito.cs.ucla.edu:27021 \
        -L 49200:pyrito.cs.ucla.edu:49200 \
        -L 49201:pyrito.cs.ucla.edu:49201 \
        -L 49202:pyrito.cs.ucla.edu:49202 \
        -L 49203:pyrito.cs.ucla.edu:49203 \
        -L 49204:pyrito.cs.ucla.edu:49204 \
        -o ServerAliveInterval=60 \
        -o ServerAliveCountMax=3 \
        -o ExitOnForwardFailure=yes \
        -o StrictHostKeyChecking=no \
        -o UserKnownHostsFile=/dev/null \
        -o LogLevel=ERROR \
        fts@scai1.cs.ucla.edu \
        > "$TUNNEL_LOG_FILE" 2>&1 &
    
    local tunnel_pid=$!
    echo $tunnel_pid > "$TUNNEL_PID_FILE"
    
    # Wait a moment and check if it started successfully
    sleep 2
    if kill -0 $tunnel_pid 2>/dev/null; then
        echo "✓ Tunnel started successfully (PID: $tunnel_pid)"
        echo "  Log file: $TUNNEL_LOG_FILE"
        return 0
    else
        echo "✗ Failed to start tunnel"
        rm -f "$TUNNEL_PID_FILE"
        return 1
    fi
}

stop_tunnel() {
    if [ -f "$TUNNEL_PID_FILE" ]; then
        local pid=$(cat "$TUNNEL_PID_FILE")
        if kill -0 $pid 2>/dev/null; then
            echo "Stopping tunnel (PID: $pid)..."
            kill $pid
            rm -f "$TUNNEL_PID_FILE"
            echo "✓ Tunnel stopped"
        else
            echo "Tunnel process not found, cleaning up PID file"
            rm -f "$TUNNEL_PID_FILE"
        fi
    else
        echo "No tunnel PID file found"
    fi
}

is_tunnel_running() {
    if [ -f "$TUNNEL_PID_FILE" ]; then
        local pid=$(cat "$TUNNEL_PID_FILE")
        if kill -0 $pid 2>/dev/null; then
            return 0
        else
            # PID file exists but process is dead
            rm -f "$TUNNEL_PID_FILE"
            return 1
        fi
    else
        return 1
    fi
}

status_tunnel() {
    if is_tunnel_running; then
        local pid=$(cat "$TUNNEL_PID_FILE")
        echo "✓ Tunnel is running (PID: $pid)"
        
        # Check if ports are actually listening
        echo "  Checking ports:"
        for port in 27020 27021; do
            if nc -z localhost $port 2>/dev/null; then
                echo "    ✓ Port $port is open"
            else
                echo "    ✗ Port $port is not responding"
            fi
        done
    else
        echo "✗ Tunnel is not running"
    fi
}

restart_tunnel() {
    stop_tunnel
    sleep 1
    start_tunnel
}

# Main command processing
case "${1:-status}" in
    start)
        start_tunnel
        ;;
    stop)
        stop_tunnel
        ;;
    restart)
        restart_tunnel
        ;;
    status)
        status_tunnel
        ;;
    *)
        echo "Usage: $0 {start|stop|restart|status}"
        echo "  start   - Start the tunnel daemon"
        echo "  stop    - Stop the tunnel daemon"
        echo "  restart - Restart the tunnel daemon"
        echo "  status  - Check tunnel status"
        exit 1
        ;;
esac

