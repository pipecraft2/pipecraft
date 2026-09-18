#!/bin/bash
# Bind-mount ownership helper for Docker and Podman (including rootless).
#
# Rootless Podman maps container UID 0 to the host user. chown $HOST_UID then
# retargets files at a subordinate UID and the host user can no longer read them.
# Skip chown in a user namespace; only retag when we are real root (Docker /
# rootful Podman).

chown_bind_to_host() {
  local target=$1
  if [ -z "${HOST_UID:-}" ] || [ -z "${HOST_GID:-}" ]; then
    return 0
  fi
  if [ ! -e "$target" ]; then
    return 0
  fi

  local mapped
  mapped=$(awk '$1==0 {print $2; exit}' /proc/self/uid_map 2>/dev/null)
  if [ -n "$mapped" ] && [ "$mapped" != "0" ]; then
    echo "Skipping chown on $target in a rootless user namespace (already mapped to the host user)"
    return 0
  fi

  echo "Setting ownership of $target to $HOST_UID:$HOST_GID"
  chown -R "$HOST_UID:$HOST_GID" "$target"
}
