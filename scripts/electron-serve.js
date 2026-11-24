"use strict";

/**
 * Cross-platform launcher for `yarn electron:serve`.
 *
 * Linux needs ELECTRON_DISABLE_SANDBOX=1 before the Electron process starts
 * (Ubuntu 24.04+ and dev builds without a setuid chrome-sandbox). The
 * `VAR=value command` form only works in sh/bash; cmd.exe on Windows treats
 * ELECTRON_DISABLE_SANDBOX as the executable name.
 */
const { spawn } = require("child_process");

const env = { ...process.env };

if (process.platform === "linux") {
  env.ELECTRON_DISABLE_SANDBOX = "1";
}

const child = spawn(
  "vue-cli-service",
  ["electron:serve", ...process.argv.slice(2)],
  {
    stdio: "inherit",
    env,
    shell: true,
  }
);

child.on("error", (err) => {
  console.error(err);
  process.exit(1);
});

child.on("exit", (code, signal) => {
  if (signal) {
    process.kill(process.pid, signal);
    return;
  }
  process.exit(code === null ? 1 : code);
});
