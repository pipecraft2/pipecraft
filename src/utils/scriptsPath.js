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

function getHostScriptsPath() {
  return path.join(getElectronApp().getPath("userData"), "service_scripts");
}

/**
 * Returns the host path to service_scripts for Docker bind mounts.
 * On Linux production (AppImage), copies bundled scripts to userData because
 * Docker cannot reliably bind-mount paths under /tmp/.mount_* FUSE mounts.
 */
function getServiceScriptsPath() {
  if (isDevelopment) {
    return slash(path.join(process.cwd(), "src/pipecraft-core/service_scripts"));
  }

  const bundled = getBundledScriptsPath();
  if (process.platform !== "linux") {
    return bundled;
  }

  const hostPath = getHostScriptsPath();
  if (!fs.existsSync(hostPath)) {
    copyDirSync(bundled, hostPath);
  }
  return hostPath;
}

module.exports = {
  getServiceScriptsPath,
  getBundledScriptsPath,
  getHostScriptsPath,
};
