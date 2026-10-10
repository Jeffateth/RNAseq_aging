#!/bin/bash
# Persistent SSH connection ("tunnel") to Euler: authenticate once, then ssh / scp / rsync reuse it for hours.
# Needs the ETH network or ETH VPN.      Usage:  EULER_USER=<nethz-username> bash euler/tunnel.sh init|start|status|stop
set -euo pipefail
: "${EULER_USER:?set EULER_USER=<your ETH username>, e.g.  EULER_USER=jyao bash euler/tunnel.sh start}"
CONF="$HOME/.ssh/euler.conf"

case "${1:-}" in
  init)    # writes ~/.ssh/euler.conf (does not touch ~/.ssh/config)
    mkdir -p "$HOME/.ssh/sockets"; chmod 700 "$HOME/.ssh"
    cat > "$CONF" <<EOC
Host euler
  HostName euler.ethz.ch
  User $EULER_USER
  ControlMaster auto
  ControlPath $HOME/.ssh/sockets/euler-%r
  ControlPersist 8h
  ServerAliveInterval 60
  ServerAliveCountMax 5
EOC
    echo "wrote $CONF"
    echo "Optional: add the line  'Include ~/.ssh/euler.conf'  at the TOP of ~/.ssh/config so plain 'ssh euler', 'scp', 'rsync' work everywhere." ;;
  start)   # opens the master connection in the background; asks for your password/key passphrase once
    [[ -f "$CONF" ]] || bash "$0" init
    ssh -F "$CONF" -fN euler && echo "tunnel up. Use:  ssh -F $CONF euler   (or 'ssh euler' after the Include line)" ;;
  status)  ssh -F "$CONF" -O check euler ;;
  stop)    ssh -F "$CONF" -O exit euler ;;
  *) echo "usage: EULER_USER=<user> $0 init|start|status|stop"; exit 1 ;;
esac
