const fs = require("fs");
const path = require("path");
const slash = require("slash");

const isDevelopment = process.env.NODE_ENV !== "production";

function getElectronApp() {
  if (process.type === "browser") {
    return require("electron").app;
  }
  return require("@electron/remote").app;
}

function copyDirSync(srcDir, destDir) {
  const entries = fs.readdirSync(srcDir, { withFileTypes: true });
  fs.mkdirSync(destDir, { recursive: true });
  for (const entry of entries) {
    const srcPath = path.join(srcDir, entry.name);
    const destPath = path.join(destDir, entry.name);
    if (entry.isDirectory()) {
      copyDirSync(srcPath, destPath);
      continue;
    }
    if (entry.isFile()) {
      fs.copyFileSync(srcPath, destPath);
    }
  }
}

function getBundledScriptsPath() {
  return path.join(
    process.resourcesPath,
    "src",
    "pipecraft-core",
    "service_scripts"
  );
}

function getUserDataPath() {
  return getElectronApp().getPath("userData");
}

/** Per-run workflow settings (NextITS, FunBarONT, OptimOTU). Not the scripts folder. */
function getWorkflowConfigPath(fileName) {
  const dir = path.join(getUserDataPath(), "workflow_configs");
  fs.mkdirSync(dir, { recursive: true });
  return path.join(dir, fileName);
}

/**
 * Params file Nextflow reads. chimera_methods has no screen control, and a
 * FASTA chimera database is always used from /dbcache/db.udb. The environment
 * variable keeps the selected file so the script can build that UDB.
 */
function applyNextItsChimeraParams(params) {
  const config = { ...params };
  const source = config.chimera_db;
  const hasDb = Boolean(source) && source !== "undefined";
  if (!hasDb) {
    delete config.chimera_db;
    config.chimera_methods = "denovo";
    return config;
  }
  config.chimera_methods = "ref,denovo";
  if (!String(source).toLowerCase().endsWith(".udb")) {
    config.chimera_db = "/dbcache/db.udb";
  }
  return config;
}

function getHostScriptsPath() {
  return path.join(getUserDataPath(), "service_scripts");
}

const VERSION_STAMP_FILE = ".pipecraft-version";

// Linux production (AppImage) uses a userData copy of service_scripts because
// the engine cannot reliably bind-mount paths under /tmp/.mount_* FUSE mounts.
function needsHostCopy() {
  return !isDevelopment && process.platform === "linux";
}

/**
 * Refresh the userData copy of service_scripts when the app version changes,
 * so updated scripts reach the containers. Call once at startup from the main
 * process, before any container runs. No-op outside Linux production.
 */
function syncServiceScripts(appVersion) {
  if (!needsHostCopy()) {
    return;
  }
  const hostPath = getHostScriptsPath();
  const stampFile = path.join(hostPath, VERSION_STAMP_FILE);
  const copiedVersion = fs.existsSync(stampFile)
    ? fs.readFileSync(stampFile, "utf8").trim()
    : "";
  if (copiedVersion === appVersion) {
    return;
  }
  fs.rmSync(hostPath, { recursive: true, force: true });
  copyDirSync(getBundledScriptsPath(), hostPath);
  // Written last so a failed copy is retried on the next launch.
  fs.writeFileSync(stampFile, appVersion);
}

/**
 * Returns the host path to service_scripts for Docker/Podman bind mounts.
 */
function getServiceScriptsPath() {
  if (isDevelopment) {
    return slash(path.join(process.cwd(), "src/pipecraft-core/service_scripts"));
  }
  if (!needsHostCopy()) {
    return getBundledScriptsPath();
  }
  const hostPath = getHostScriptsPath();
  if (!fs.existsSync(hostPath)) {
    copyDirSync(getBundledScriptsPath(), hostPath);
  }
  return hostPath;
}

module.exports = {
  syncServiceScripts,
  getServiceScriptsPath,
  getBundledScriptsPath,
  getHostScriptsPath,
  getUserDataPath,
  getWorkflowConfigPath,
  applyNextItsChimeraParams,
};
