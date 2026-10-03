#!/usr/bin/env python3
"""SSH helper for SDSC Expanse (account cla361, user fsun3) for the ML-LUCJ project.

Expanse requires TOTP 2FA. This helper authenticates once and keeps an SSH
ControlMaster socket open (default 8 h), so later commands, rsyncs and status
checks reuse the authenticated connection without another TOTP code.

Credentials are NOT stored here: the TOTP seed is read from
EXPANSE_TOTP_SECRET / EXPANSE_TOTP_SECRET_FILE (default
~/.config/fermi_arc/expanse_totp_secret, mode 0600, shared with the user's other
project). The last-used TOTP window is recorded in the sibling
expanse_last_totp_counter file under an exclusive flock, the same protocol the
other project's helper uses, so concurrent helpers never reuse a code.

Usage (from scai1/scai2):
    python3 expanse.py connect              # open (or reuse) the master connection
    python3 expanse.py check                # is the master alive?
    python3 expanse.py run "squeue -u fsun3"
    python3 expanse.py put LOCAL REMOTE     # rsync to Expanse (dirs: add trailing /)
    python3 expanse.py get REMOTE LOCAL     # rsync from Expanse
    python3 expanse.py status               # queue + allocation balance
    python3 expanse.py close                # close the master connection
"""
from __future__ import annotations

import fcntl
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

HOST = os.environ.get("EXPANSE_SSH_HOST", "fsun3@login.expanse.sdsc.edu")
KEY = Path(os.environ.get("EXPANSE_SSH_KEY", "~/.ssh/id_rsa")).expanduser()
SECRET_FILE = Path(os.environ.get(
    "EXPANSE_TOTP_SECRET_FILE", "~/.config/fermi_arc/expanse_totp_secret")).expanduser()
COUNTER_FILE = SECRET_FILE.with_name("expanse_last_totp_counter")
CTL_DIR = Path("~/.ssh/cm").expanduser()
CTL = CTL_DIR / "expanse-fsun3"
PERSIST = os.environ.get("EXPANSE_CONTROL_PERSIST", "8h")
ACCOUNT = "cla361"
REMOTE_ROOT = "/expanse/lustre/scratch/fsun3/temp_project/ml_lucj"

BASE_OPTS = ["-i", str(KEY), "-o", "StrictHostKeyChecking=accept-new",
             "-o", "ServerAliveInterval=60", "-o", "ServerAliveCountMax=5"]
MUX_OPTS = ["-o", f"ControlPath={CTL}", "-o", "ControlMaster=no", "-o", "BatchMode=yes"]


def _secret() -> str:
    s = os.environ.get("EXPANSE_TOTP_SECRET", "").strip()
    if not s:
        try:
            s = SECRET_FILE.read_text().strip()
        except FileNotFoundError as exc:
            raise SystemExit(f"TOTP seed not found ({SECRET_FILE}); see the Expanse guide.") from exc
    if not s:
        raise SystemExit("TOTP seed is empty")
    return s


def _reserve_totp(min_remaining: float = 12.0) -> str:
    """A code from a TOTP window no helper has used yet (shared flock-protected counter)."""
    import pyotp
    totp = pyotp.TOTP(_secret())
    COUNTER_FILE.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(COUNTER_FILE, os.O_RDWR | os.O_CREAT, 0o600)
    with os.fdopen(fd, "r+") as f:
        os.chmod(COUNTER_FILE, 0o600)
        fcntl.flock(f, fcntl.LOCK_EX)
        while True:
            now = time.time()
            counter = int(now // totp.interval)
            remaining = totp.interval - now % totp.interval
            f.seek(0)
            raw = f.read().strip()
            last = int(raw) if raw else -1
            if counter > last and remaining >= min_remaining:
                break
            time.sleep(remaining + 0.5)
        f.seek(0)
        f.truncate()
        f.write(f"{counter}\n")
        f.flush()
        os.fsync(f.fileno())
        return totp.at(int(time.time()))


def master_alive() -> bool:
    if not CTL.exists():
        return False
    r = subprocess.run(["ssh", "-o", f"ControlPath={CTL}", "-O", "check", HOST],
                       capture_output=True, text=True)
    return r.returncode == 0


def connect(timeout: int = 120) -> None:
    """Open a background ControlMaster connection (one TOTP code)."""
    if master_alive():
        print("master connection already alive")
        return
    import pexpect
    CTL_DIR.mkdir(mode=0o700, parents=True, exist_ok=True)
    os.chmod(CTL_DIR, 0o700)
    if CTL.exists():
        CTL.unlink()  # stale socket
    # ControlMaster=auto + ControlPersist: the first session authenticates, runs a
    # marker command, and OpenSSH forks the master into the background (setsid),
    # where it survives this pexpect pty being closed.
    cmd = ["ssh", *BASE_OPTS, "-o", "PreferredAuthentications=publickey,keyboard-interactive",
           "-o", f"ControlPath={CTL}", "-o", "ControlMaster=auto",
           "-o", f"ControlPersist={PERSIST}", HOST, "echo EXPANSE_MASTER_OK"]
    child = pexpect.spawn(cmd[0], cmd[1:], timeout=timeout, encoding="utf-8")
    sent = []
    while True:
        i = child.expect([r"TOTP code for \w+:", r"\(no-PIN\) Yubi for \w+:", r"Password:",
                          r"Permission denied", r"EXPANSE_MASTER_OK", pexpect.EOF, pexpect.TIMEOUT])
        if i == 0:
            code = _reserve_totp()
            sent.append(code)
            child.sendline(code)
        elif i in (1, 2):
            child.sendline("")
        elif i == 3:
            raise SystemExit("authentication failed (Permission denied)")
        elif i in (4, 5):
            break
        else:
            raise SystemExit("timed out while authenticating")
    child.close(force=True)
    for _ in range(20):
        if master_alive():
            print(f"master connection up (ControlPersist={PERSIST})")
            return
        time.sleep(0.5)
    out = child.before or ""
    for c in sent:
        out = out.replace(c, "[TOTP]")
    raise SystemExit(f"master connection did not come up: {out[-300:]}")


def _ensure() -> None:
    if not master_alive():
        connect()


def run(command: str, check: bool = False, capture: bool = False, timeout: int | None = None):
    _ensure()
    r = subprocess.run(["ssh", *BASE_OPTS, *MUX_OPTS, HOST, command],
                       capture_output=capture, text=True, timeout=timeout)
    if check and r.returncode != 0:
        raise SystemExit(f"remote command failed ({r.returncode}): {command}\n{(r.stderr or '')[-500:]}")
    return r


def _rsync(src: str, dst: str, extra: list[str] | None = None) -> int:
    _ensure()
    ssh = " ".join(shlex.quote(x) for x in ["ssh", *BASE_OPTS, *MUX_OPTS])
    cmd = ["rsync", "-az", "--partial", "-e", ssh, *(extra or []), src, dst]
    return subprocess.run(cmd).returncode


def put(local: str, remote: str, extra: list[str] | None = None) -> int:
    return _rsync(local, f"{HOST}:{remote}", extra)


def get(remote: str, local: str, extra: list[str] | None = None) -> int:
    return _rsync(f"{HOST}:{remote}", local, extra)


def status() -> None:
    run(f"squeue -u fsun3 -o '%.12i %.9P %.24j %.8T %.10M %.10l %.5C %R' ; "
        f"echo; expanse-client user -r expanse -p 2>/dev/null | head -20")


def close() -> None:
    if CTL.exists():
        subprocess.run(["ssh", "-o", f"ControlPath={CTL}", "-O", "exit", HOST])
    print("closed")


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 1
    cmd, rest = argv[0], argv[1:]
    if cmd == "connect":
        connect()
    elif cmd == "check":
        print("alive" if master_alive() else "down")
    elif cmd == "run":
        return run(" ".join(rest)).returncode
    elif cmd == "put":
        return put(rest[0], rest[1], rest[2:])
    elif cmd == "get":
        return get(rest[0], rest[1], rest[2:])
    elif cmd == "status":
        status()
    elif cmd == "close":
        close()
    else:
        print(__doc__)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
