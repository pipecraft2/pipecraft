const { execFileSync } = require("child_process");
const fs = require("fs");
const os = require("os");
const path = require("path");

function parseDockerHostToOptions(dockerHost) {
  if (!dockerHost) {
    throw new Error("Docker context returned an empty endpoint.");
  }

  if (dockerHost.startsWith("unix://")) {
    return { socketPath: dockerHost.replace("unix://", "") };
  }

  if (dockerHost.startsWith("npipe://")) {
    // Docker CLI reports Windows named pipes as npipe://..., but Node expects \\.\pipe\<name>.
    // Example: npipe:////./pipe/dockerDesktopLinuxEngine -> \\.\pipe\dockerDesktopLinuxEngine
    const pipeMarker = "/pipe/";
    const idx = dockerHost.indexOf(pipeMarker);
    if (idx === -1) {
      throw new Error(`Unsupported npipe endpoint from context: ${dockerHost}`);
    }
    const pipeName = dockerHost.slice(idx + pipeMarker.length);
    return { socketPath: `\\\\.\\pipe\\${pipeName}` };
  }

  if (dockerHost.startsWith("tcp://")) {
    const endpoint = new URL(dockerHost);
    return {
      host: endpoint.hostname,
      port: Number(endpoint.port || 2375),
      protocol: "http",
    };
  }

  throw new Error(`Unsupported Docker endpoint from context: ${dockerHost}`);
}

function resolveExistingSocketPath(preferredSocketPath) {
  const candidates = [];

  if (preferredSocketPath) {
    candidates.push(preferredSocketPath);
  }

  // Common locations across Linux + Docker Desktop for macOS.
  candidates.push("/var/run/docker.sock");
  candidates.push(path.join(os.homedir(), ".docker", "run", "docker.sock"));
  candidates.push(path.join(os.homedir(), ".docker", "desktop", "docker.sock"));

  const unique = [...new Set(candidates)].filter(Boolean);
  const existing = unique.find((p) => {
    try {
      return fs.existsSync(p);
    } catch {
      return false;
    }
  });

  return existing || preferredSocketPath;
}

function readContextValue(args) {
  return execFileSync("docker", args, {
    encoding: "utf8",
    windowsHide: true,
    stdio: ["ignore", "pipe", "pipe"],
  }).trim();
}

function getDockerodeOptionsFromContextSync() {
  try {
    const contextName = readContextValue(["context", "show"]);
    if (!contextName) {
      throw new Error("No active Docker context found.");
    }

    const dockerHost = readContextValue([
      "context",
      "inspect",
      contextName,
      "--format",
      '{{ (index .Endpoints "docker").Host }}',
    ]);

    const opts = parseDockerHostToOptions(dockerHost);
    if (opts?.socketPath) {
      return { ...opts, socketPath: resolveExistingSocketPath(opts.socketPath) };
    }
    return opts;
  } catch (error) {
    if (error?.code === "ENOENT") {
      // In packaged Electron apps, PATH may not include the docker CLI even if Docker Desktop is running.
      // Fall back to probing common socket paths so we can still connect via Dockerode.
      const fallbackSocketPath = resolveExistingSocketPath();
      if (fallbackSocketPath) {
        return { socketPath: fallbackSocketPath };
      }
      throw new Error("Docker CLI is not installed or not available in PATH.");
    }

    const stderr = error?.stderr ? String(error.stderr).trim() : "";
    const extra = stderr ? ` (${stderr})` : "";
    throw new Error(`Failed to resolve Docker endpoint from Docker context${extra}`);
  }
}

module.exports = { getDockerodeOptionsFromContextSync };
