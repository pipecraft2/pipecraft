#!/bin/bash
# Reusable database index cache (BLAST databases, vsearch .udb files).
#
# PipeCraft mounts one writable cache folder per database at /dbcache
# (host: <PipeCraft app data>/db_cache/<format>/<entry>, see src/utils/dbCache.js).
# An index built there is reused by later runs, in any project, as long as the
# source database file is unchanged. The ".complete" marker is written only
# after a successful build, so an interrupted build is redone next time.
#
# Usage:
#   cache=$(db_cache_dir)
#   if ! db_cache_is_complete "$cache"; then
#       db_cache_reset "$cache"
#       <build index into $cache> || exit 1
#       db_cache_mark_complete "$cache"
#   fi

DB_CACHE_MOUNT="/dbcache"

# Prints the cache folder; falls back to a throwaway temp folder when PipeCraft
# did not mount one (the index is then rebuilt on every run).
db_cache_dir() {
    if [[ -d "$DB_CACHE_MOUNT" ]] && [[ -w "$DB_CACHE_MOUNT" ]]; then
        echo "$DB_CACHE_MOUNT"
    else
        mktemp -d
    fi
}

db_cache_is_complete() {
    [[ -f "$1/.complete" ]]
}

# Removes leftovers of an interrupted build; keeps PipeCraft's info.json.
db_cache_reset() {
    find "$1" -mindepth 1 -maxdepth 1 ! -name info.json -exec rm -rf {} +
}

db_cache_mark_complete() {
    touch "$1/.complete"
}
