/**
 * Database index cache shared by all projects (scripts side:
 * service_scripts/submodules/db_cache.sh).
 *
 * Layout: <userData>/db_cache/<format>/<entry>/
 *   index files   - written by the container (BLAST database, vsearch .udb)
 *   .complete     - written by the container after a successful build
 *   info.json     - written here; source file details for the cache manager
 *
 * An entry is identified by the source file's name, size, modification time
 * and the image that builds the index, so an edited database or a new tool
 * version gets a fresh index instead of a stale one.
 */
const crypto = require("crypto");
const fs = require("fs");
const path = require("path");
const slash = require("slash");
const { getUserDataPath } = require("./scriptsPath");
const { isFileInputInUse } = require("./fileMounts");

const CACHE_MOUNT = "/dbcache";
const INFO_FILE = "info.json";
const COMPLETE_MARKER = ".complete";

const CACHED_INPUT_BY_SCRIPT = {
  "taxonomy_BLAST.sh": { input: "database_file", format: "blast" },
  "taxonomy_sintax.sh": { input: "database", format: "udb" },
  "NextITS_Pipeline.sh": { input: "chimera_db", format: "udb" },
};

const FORMAT_LABELS = { blast: "BLAST", udb: "UDB" };

// Already an index; scripts use these files directly.
const PREBUILT_EXTENSIONS = new Set([
  ".udb", ".ndb", ".nhr", ".nin", ".nsq", ".not", ".ntf", ".nto", ".nal",
]);

function getCacheRoot() {
  return path.join(getUserDataPath(), "db_cache");
}

function entryName(hostFile, stat, imageName) {
  const base = path
    .basename(hostFile)
    .replace(/[^A-Za-z0-9._-]/g, "_")
    .slice(0, 60);
  const hash = crypto
    .createHash("sha1")
    .update(`${stat.size}:${Math.floor(stat.mtimeMs)}:${imageName || ""}`)
    .digest("hex")
    .slice(0, 10);
  return `${base}_${hash}`;
}

function readInfo(entryDir) {
  try {
    return JSON.parse(fs.readFileSync(path.join(entryDir, INFO_FILE), "utf8"));
  } catch {
    return {};
  }
}

// BLAST databases built next to the fasta by older PipeCraft versions are used in place.
function hasLegacyBlastIndex(hostFile) {
  return [".nhr", ".nal"].some((ext) => fs.existsSync(`${hostFile}${ext}`));
}

/**
 * Prepares the cache entry for a service's selected database and returns the
 * bind that mounts it at /dbcache, or null when nothing needs caching.
 */
function prepareDbCacheBind(element, inputs, isHidden) {
  const cached = CACHED_INPUT_BY_SCRIPT[element.scriptName];
  if (!cached) {
    return null;
  }
  const input = inputs.find(
    (candidate) => candidate.name === cached.input && isFileInputInUse(candidate, isHidden)
  );
  if (!input) {
    return null;
  }
  const hostFile = slash(String(input.value));
  if (
    PREBUILT_EXTENSIONS.has(path.extname(hostFile).toLowerCase()) ||
    !fs.existsSync(hostFile) ||
    (cached.format === "blast" && hasLegacyBlastIndex(hostFile))
  ) {
    return null;
  }

  const stat = fs.statSync(hostFile);
  const entryDir = path.join(
    getCacheRoot(),
    cached.format,
    entryName(hostFile, stat, element.imageName)
  );
  fs.mkdirSync(entryDir, { recursive: true });

  const now = new Date().toISOString();
  const info = {
    format: cached.format,
    sourceName: path.basename(hostFile),
    sourcePath: hostFile,
    sourceSize: stat.size,
    sourceMtimeMs: Math.floor(stat.mtimeMs),
    image: element.imageName || "",
    createdAt: readInfo(entryDir).createdAt || now,
    lastUsedAt: now,
  };
  fs.writeFileSync(path.join(entryDir, INFO_FILE), JSON.stringify(info, null, 2));
  return `${slash(entryDir)}:${CACHE_MOUNT}`;
}

function dirSize(dir) {
  return fs.readdirSync(dir, { withFileTypes: true }).reduce((total, entry) => {
    const fullPath = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      return total + dirSize(fullPath);
    }
    return total + fs.statSync(fullPath).size;
  }, 0);
}

// "ok" | "incomplete" | "missing" (source deleted/moved) | "changed" (source edited)
function entryStatus(entryDir, info) {
  if (!fs.existsSync(path.join(entryDir, COMPLETE_MARKER))) {
    return "incomplete";
  }
  if (!info.sourcePath || !fs.existsSync(info.sourcePath)) {
    return "missing";
  }
  const stat = fs.statSync(info.sourcePath);
  if (stat.size !== info.sourceSize || Math.floor(stat.mtimeMs) !== info.sourceMtimeMs) {
    return "changed";
  }
  return "ok";
}

/** All cache entries, most recently used first. */
function listCacheEntries() {
  const root = getCacheRoot();
  if (!fs.existsSync(root)) {
    return [];
  }
  const entries = [];
  fs.readdirSync(root, { withFileTypes: true })
    .filter((formatDir) => formatDir.isDirectory())
    .forEach((formatDir) => {
      const formatPath = path.join(root, formatDir.name);
      fs.readdirSync(formatPath, { withFileTypes: true })
        .filter((entry) => entry.isDirectory())
        .forEach((entry) => {
          const entryDir = path.join(formatPath, entry.name);
          const info = readInfo(entryDir);
          entries.push({
            id: entryDir,
            name: info.sourceName || entry.name,
            format: FORMAT_LABELS[formatDir.name] || formatDir.name,
            sizeBytes: dirSize(entryDir),
            lastUsedAt: info.lastUsedAt || "",
            sourcePath: info.sourcePath || "",
            status: entryStatus(entryDir, info),
          });
        });
    });
  return entries.sort((a, b) => b.lastUsedAt.localeCompare(a.lastUsedAt));
}

/** Deletes the given entries (ids from listCacheEntries); ignores paths outside the cache. */
function removeCacheEntries(ids) {
  const root = path.resolve(getCacheRoot());
  ids
    .map((id) => path.resolve(id))
    .filter((dir) => dir.startsWith(root + path.sep))
    .forEach((dir) => fs.rmSync(dir, { recursive: true, force: true }));
}

module.exports = {
  getCacheRoot,
  listCacheEntries,
  prepareDbCacheBind,
  removeCacheEntries,
};
