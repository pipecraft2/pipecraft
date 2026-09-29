/**
 * chrome-sandbox cannot be root+setuid on a FUSE mount, and Ubuntu 24.04+
 * often blocks the userns fallback. Wrap the Electron binary so every launch
 * sets ELECTRON_DISABLE_SANDBOX and passes --no-sandbox, then remove
 * chrome-sandbox.
 */
const fs = require("fs");
const path = require("path");

exports.default = async function afterPack(context) {
  if (context.electronPlatformName !== "linux") {
    return;
  }

  const { appOutDir, packager } = context;
  const exeName = packager.executableName;
  const exePath = path.join(appOutDir, exeName);
  const wrappedPath = path.join(appOutDir, `${exeName}.bin`);
  const chromeSandbox = path.join(appOutDir, "chrome-sandbox");

  if (fs.existsSync(chromeSandbox)) {
    fs.unlinkSync(chromeSandbox);
  }

  if (!fs.existsSync(exePath)) {
    console.warn(`[afterPack] Linux executable not found: ${exePath}`);
    return;
  }

  // Already wrapped from a previous pack in the same out dir
  if (fs.existsSync(wrappedPath)) {
    return;
  }

  fs.renameSync(exePath, wrappedPath);
  const wrapper = `#!/bin/bash
DIR="$(cd "$(dirname "$0")" && pwd)"
export ELECTRON_DISABLE_SANDBOX=1
exec "$DIR/${exeName}.bin" --no-sandbox "$@"
`;
  fs.writeFileSync(exePath, wrapper, { mode: 0o755 });
  console.log(`[afterPack] Wrapped ${exeName} with --no-sandbox launcher`);
};
