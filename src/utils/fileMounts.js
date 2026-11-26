/**
 * Mounts for user-selected files (databases, tables, fasta files, ...).
 *
 * Each selected file is mounted at /extraFiles/<input name>/<file name>, and
 * the script receives that container path in its env variable, so scripts
 * never need to know mount slots or rebuild paths.
 *
 * Mount modes (default: the single file, read-only):
 *   "folder"    - the file's folder, read-only (tools that need sibling files)
 *   "folder-rw" - the file's folder, writable (tools that write indexes next to the file)
 * Modes are keyed by script name so they also apply to loaded older configs.
 */
const fs = require("fs");
const path = require("path");
const slash = require("slash");

const MOUNT_ROOT = "/extraFiles";

const MOUNT_MODES_BY_SCRIPT = {
  "taxonomy_BLAST.sh": { database_file: "folder" },
  "taxonomy_RDP.sh": { database: "folder" },
  "blasch.sh": { reference_db: "folder" },
};

function hasFileValue(input) {
  return Boolean(input.value) && input.value !== "undefined";
}

function isFileInputInUse(input, isHidden) {
  if (input.type === "boolfile") {
    return input.active === true && hasFileValue(input);
  }
  if (input.type !== "file" || !hasFileValue(input)) {
    return false;
  }
  return input.depends_on === undefined || !isHidden(input);
}

/**
 * Selected file inputs with their host and container paths.
 * @param {Array} inputs - service Inputs + extraInputs
 * @param {Function} isHidden - store getter check_depends_on
 * @returns {Array<{input, hostFile, mountDir, containerPath}>}
 */
function collectFileInputs(inputs, isHidden) {
  const usedDirs = new Set();
  return inputs
    .filter((input) => isFileInputInUse(input, isHidden))
    .map((input) => {
      let mountDir = `${MOUNT_ROOT}/${input.name}`;
      for (let n = 2; usedDirs.has(mountDir); n++) {
        mountDir = `${MOUNT_ROOT}/${input.name}_${n}`;
      }
      usedDirs.add(mountDir);
      const hostFile = slash(String(input.value));
      return {
        input,
        hostFile,
        mountDir,
        containerPath: `${mountDir}/${path.posix.basename(hostFile)}`,
      };
    });
}

/** Map of input object -> container path, for building env variables. */
function containerPathsByInput(inputs, isHidden) {
  return new Map(
    collectFileInputs(inputs, isHidden).map((entry) => [entry.input, entry.containerPath])
  );
}

function bindFor(entry, mode) {
  const hostDir = path.posix.dirname(entry.hostFile);
  if (mode === "folder-rw") {
    return `${hostDir}:${entry.mountDir}`;
  }
  if (mode === "folder") {
    return `${hostDir}:${entry.mountDir}:ro`;
  }
  return `${entry.hostFile}:${entry.containerPath}:ro`;
}

/**
 * Bind specs for the selected files of one service.
 * Throws when a selected file no longer exists: Docker would otherwise
 * create an empty folder at that path.
 */
function buildFileBinds(inputs, isHidden, scriptName) {
  const modes = MOUNT_MODES_BY_SCRIPT[scriptName] || {};
  const entries = collectFileInputs(inputs, isHidden);

  const missing = entries.filter((entry) => !fs.existsSync(entry.hostFile));
  if (missing.length) {
    const list = missing
      .map((entry) => `${entry.input.name.replace(/_/g, " ")}: ${entry.hostFile}`)
      .join("\n");
    throw new Error(`Selected file not found. Select it again and restart.\n${list}`);
  }

  return entries.map((entry) => bindFor(entry, modes[entry.input.name]));
}

module.exports = {
  buildFileBinds,
  containerPathsByInput,
  isFileInputInUse,
};
