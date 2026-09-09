/**
 * Docker + Podman helpers for PipeCraft.
 *
 * dockerode talks to both engines over the same Docker-compatible API.
 * This module finds the right socket/pipe, remembers the choice, and
 * fixes Windows Podman bind mounts before a container starts.
 */
const { execFile, execFileSync, spawn } = require("child_process");
const fs = require("fs");
const os = require("os");
const path = require("path");
const { promisify } = require("util");

const execFileAsync = promisify(execFile);

// Saved in localStorage: "auto" | "docker" | "podman"
const PREFERENCE_KEY = "pipecraft.containerRuntime";
const PREFERENCES = ["auto", "docker", "podman"];
const CLI_TIMEOUT_MS = 8000;
const MACHINE_TIMEOUT_MS = 120000;

// Last engine we successfully connected to (socket options, engine name, etc.)
let cachedRuntime = null;
// Background `podman system service` process, if we had to start one ourselves.
let podmanServiceChild = null;

function readRuntimePreference() {
  try {
    if (typeof localStorage !== "undefined") {
      const stored = localStorage.getItem(PREFERENCE_KEY);
      if (PREFERENCES.includes(stored)) {
        return stored;
      }
    }
  } catch {
    // localStorage is unavailable in the Electron main process.
  }
  return "auto";
}

function persistRuntimePreference(preference) {
  if (!PREFERENCES.includes(preference)) {
    throw new Error(`Unsupported container runtime preference: ${preference}`);
  }
  try {
    if (typeof localStorage !== "undefined") {
      localStorage.setItem(PREFERENCE_KEY, preference);
    }
  } catch {
    // Ignore persistence failures; in-memory preference still applies.
  }
}

// Extra folders to search when PATH does not include Docker/Podman.
function extraBinDirs() {
  const home = os.homedir();
  const programFiles = process.env.ProgramFiles || "C:\\Program Files";
  const localAppData = process.env.LOCALAPPDATA || path.join(home, "AppData", "Local");

  if (process.platform === "win32") {
    return [
      path.join(programFiles, "Docker", "Docker", "resources", "bin"),
      path.join(programFiles, "RedHat", "Podman"),
      path.join(programFiles, "Podman"),
      path.join(localAppData, "Programs", "Podman"),
      path.join(localAppData, "podman"),
      path.join(programFiles, "Podman Desktop", "resources", "podman"),
      path.join(localAppData, "Programs", "podman-desktop", "resources", "podman"),
    ];
  }

  if (process.platform === "darwin") {
    return [
      "/opt/homebrew/bin",
      "/usr/local/bin",
      "/opt/podman/bin",
      "/usr/local/podman/bin",
      path.join(home, ".local", "bin"),
      "/Applications/Docker.app/Contents/Resources/bin",
      "/Applications/Podman Desktop.app/Contents/Resources/podman",
    ];
  }

  return ["/usr/bin", "/usr/local/bin", "/opt/podman/bin", path.join(home, ".local", "bin")];
}

// Find `docker` / `podman` on PATH, then fall back to known install folders.
function findExecutable(commandName) {
  const exe =
    process.platform === "win32" && !commandName.endsWith(".exe")
      ? `${commandName}.exe`
      : commandName;
  const whichCmd = process.platform === "win32" ? "where" : "which";

  try {
    const result = execFileSync(whichCmd, [commandName], {
      encoding: "utf8",
      windowsHide: true,
      timeout: 3000,
      stdio: ["ignore", "pipe", "pipe"],
    })
      .trim()
      .split(/\r?\n/)[0]
      .trim();
    if (result && fs.existsSync(result)) {
      return result;
    }
  } catch {
    // Fall through to known install locations. Packaged Electron apps often
    // have a stripped PATH that does not include Docker/Podman.
  }

  for (const dir of extraBinDirs()) {
    const candidate = path.join(dir, exe);
    if (fs.existsSync(candidate)) {
      return candidate;
    }
  }

  return null;
}

function runCliSync(binaryPath, args, timeout = CLI_TIMEOUT_MS) {
  return execFileSync(binaryPath, args, {
    encoding: "utf8",
    windowsHide: true,
    timeout,
    stdio: ["ignore", "pipe", "pipe"],
  }).trim();
}

async function runCli(binaryPath, args, timeout = CLI_TIMEOUT_MS) {
  const { stdout } = await execFileAsync(binaryPath, args, {
    encoding: "utf8",
    windowsHide: true,
    timeout,
  });
  return String(stdout || "").trim();
}

// Turn DOCKER_HOST / CONTAINER_HOST style URLs into dockerode options.
function parseHostToOptions(dockerHost) {
  if (!dockerHost) {
    throw new Error("Container engine returned an empty endpoint.");
  }

  if (dockerHost.startsWith("unix://")) {
    return { socketPath: dockerHost.replace("unix://", "") };
  }

  if (dockerHost.startsWith("npipe://")) {
    // Docker/Podman CLIs report Windows named pipes as npipe://..., but Node
    // expects \\.\pipe\<name>.
    // Example: npipe:////./pipe/dockerDesktopLinuxEngine -> \\.\pipe\dockerDesktopLinuxEngine
    const pipeMarker = "/pipe/";
    const idx = dockerHost.indexOf(pipeMarker);
    if (idx === -1) {
      throw new Error(`Unsupported npipe endpoint: ${dockerHost}`);
    }
    const pipeName = dockerHost.slice(idx + pipeMarker.length);
    return { socketPath: `\\\\.\\pipe\\${pipeName}` };
  }

  if (dockerHost.startsWith("tcp://") || dockerHost.startsWith("http://")) {
    const endpoint = new URL(dockerHost);
    return {
      host: endpoint.hostname,
      port: Number(endpoint.port || 2375),
      protocol: "http",
    };
  }

  if (dockerHost.startsWith("https://")) {
    const endpoint = new URL(dockerHost);
    return {
      host: endpoint.hostname,
      port: Number(endpoint.port || 2376),
      protocol: "https",
    };
  }

  if (dockerHost.startsWith("ssh://")) {
    throw new Error("SSH container endpoints are not supported. Use a local socket or named pipe.");
  }

  if (dockerHost.startsWith("\\\\.\\pipe\\") || dockerHost.startsWith("//./pipe/")) {
    return { socketPath: dockerHost.replace("//./pipe/", "\\\\.\\pipe\\") };
  }

  if (dockerHost.startsWith("/") || dockerHost.endsWith(".sock")) {
    return { socketPath: dockerHost };
  }

  throw new Error(`Unsupported container endpoint: ${dockerHost}`);
}

function pathLooksUsable(candidate) {
  if (!candidate) {
    return false;
  }
  try {
    return fs.existsSync(candidate);
  } catch {
    return false;
  }
}

function currentUid() {
  try {
    if (typeof process.getuid === "function") {
      return process.getuid();
    }
  } catch {
    return null;
  }
  return null;
}

function isRootUser() {
  return currentUid() === 0;
}

// Rootless Podman socket for the current Linux user (not /run/podman/...).
function defaultRootlessPodmanSocket() {
  const xdg = process.env.XDG_RUNTIME_DIR;
  if (xdg) {
    return path.join(xdg, "podman", "podman.sock");
  }
  const uid = currentUid();
  if (uid != null) {
    return path.join("/run/user", String(uid), "podman", "podman.sock");
  }
  return null;
}

// System-wide Podman socket (only useful when running as root).
function defaultRootfulPodmanSocket() {
  return "/run/podman/podman.sock";
}

function preferredLocalPodmanSocket(rootless) {
  if (isRootUser() || rootless === false) {
    return defaultRootfulPodmanSocket();
  }
  return defaultRootlessPodmanSocket() || defaultRootfulPodmanSocket();
}

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

async function waitForPath(targetPath, timeoutMs) {
  if (!targetPath) {
    return false;
  }
  const started = Date.now();
  while (Date.now() - started < timeoutMs) {
    if (pathLooksUsable(targetPath)) {
      return true;
    }
    await sleep(200);
  }
  return pathLooksUsable(targetPath);
}

function collectSocketsUnder(dir, fileName = "podman.sock") {
  const found = [];
  if (!dir || !fs.existsSync(dir)) {
    return found;
  }

  try {
    const sockAtRoot = path.join(dir, fileName);
    if (fs.existsSync(sockAtRoot)) {
      found.push(sockAtRoot);
    }

    const entries = fs.readdirSync(dir, { withFileTypes: true });
    for (const entry of entries) {
      if (!entry.isDirectory()) {
        continue;
      }
      const nested = path.join(dir, entry.name, fileName);
      if (fs.existsSync(nested)) {
        found.push(nested);
      }
    }
  } catch {
    return found;
  }

  return found;
}

function uniquePaths(candidates) {
  return [...new Set(candidates.filter(Boolean))];
}

function firstUsablePath(candidates, fallback) {
  const existing = uniquePaths(candidates).find((candidate) => pathLooksUsable(candidate));
  return existing || fallback || null;
}

// Possible Docker sockets / Windows named pipes to try, in priority order.
function dockerSocketCandidates(preferredSocketPath) {
  const home = os.homedir();
  const candidates = [];

  if (preferredSocketPath) {
    candidates.push(preferredSocketPath);
  }

  candidates.push("/var/run/docker.sock");
  candidates.push("/run/docker.sock");
  candidates.push(path.join(home, ".docker", "run", "docker.sock"));
  candidates.push(path.join(home, ".docker", "desktop", "docker.sock"));
  candidates.push(path.join(home, ".colima", "default", "docker.sock"));
  candidates.push(path.join(home, ".colima", "docker.sock"));

  if (process.platform === "win32") {
    candidates.push("\\\\.\\pipe\\dockerDesktopLinuxEngine");
    candidates.push("\\\\.\\pipe\\docker_engine");
    candidates.push("\\\\.\\pipe\\docker_engine_linux");
  }

  return uniquePaths(candidates);
}

// Possible Podman sockets / pipes. Prefer rootless; only probe rootful as root.
function podmanSocketCandidates(preferredSocketPath) {
  const home = os.homedir();
  const candidates = [];

  if (preferredSocketPath) {
    candidates.push(preferredSocketPath);
  }

  const rootlessSocket = defaultRootlessPodmanSocket();
  if (rootlessSocket) {
    candidates.push(rootlessSocket);
  }

  // Rootful socket is not the default for a normal user session.
  if (isRootUser()) {
    candidates.push(defaultRootfulPodmanSocket());
    candidates.push("/var/run/podman/podman.sock");
  }

  candidates.push(path.join(home, ".local", "share", "containers", "podman", "machine", "podman.sock"));
  candidates.push(path.join(home, ".config", "containers", "podman", "machine", "podman.sock"));

  const machineRoots = [
    path.join(home, ".local", "share", "containers", "podman", "machine"),
    path.join(home, ".config", "containers", "podman", "machine"),
    path.join(home, "Library", "Application Support", "containers", "podman", "machine"),
  ];
  for (const root of machineRoots) {
    candidates.push(...collectSocketsUnder(root));
    candidates.push(...collectSocketsUnder(path.join(root, "qemu")));
    candidates.push(...collectSocketsUnder(path.join(root, "applehv")));
    candidates.push(...collectSocketsUnder(path.join(root, "libkrun")));
    candidates.push(...collectSocketsUnder(path.join(root, "wsl")));
    candidates.push(...collectSocketsUnder(path.join(root, "hyperv")));
  }

  if (process.platform === "win32") {
    candidates.push("\\\\.\\pipe\\podman-machine-default");
    candidates.push("\\\\.\\pipe\\podman");
  }

  return uniquePaths(candidates);
}

// Honor env overrides: Podman prefers CONTAINER_HOST, Docker uses DOCKER_HOST.
function optionsFromEnv(engine) {
  const envHost =
    engine === "podman"
      ? process.env.CONTAINER_HOST || process.env.DOCKER_HOST
      : process.env.DOCKER_HOST;
  if (!envHost) {
    return null;
  }
  try {
    return parseHostToOptions(envHost);
  } catch {
    return null;
  }
}

// Ask the Docker CLI which context/socket is currently active.
function resolveDockerOptionsFromCli(dockerBin) {
  const contextName = runCliSync(dockerBin, ["context", "show"]);
  if (!contextName) {
    throw new Error("No active Docker context found.");
  }

  const dockerHost = runCliSync(dockerBin, [
    "context",
    "inspect",
    contextName,
    "--format",
    '{{ (index .Endpoints "docker").Host }}',
  ]);

  return parseHostToOptions(dockerHost);
}

function parseJsonSafe(raw) {
  if (!raw) {
    return null;
  }
  try {
    return JSON.parse(raw);
  } catch {
    return null;
  }
}

function isRootlessFromPodmanInfo(info) {
  if (!info || typeof info !== "object") {
    return false;
  }
  if (info.host && typeof info.host.security === "object") {
    return Boolean(info.host.security.rootless);
  }
  return false;
}

function readPodmanInfo(podmanBin) {
  try {
    return parseJsonSafe(runCliSync(podmanBin, ["info", "--format", "json"]));
  } catch {
    return null;
  }
}

// Read the API socket path from `podman info` (most reliable source).
function remoteSocketFromInfo(info) {
  const remote = info?.host?.remoteSocket;
  if (!remote || !remote.path) {
    return { path: null, exists: false };
  }

  let socketPath = null;
  try {
    socketPath = parseHostToOptions(remote.path).socketPath || null;
  } catch {
    socketPath = remote.path.startsWith("/") ? remote.path : null;
  }

  if (!socketPath) {
    return { path: null, exists: false };
  }

  if (process.platform === "win32" && socketPath.includes("\\pipe\\")) {
    return { path: socketPath, exists: remote.exists !== false };
  }

  return {
    path: socketPath,
    exists: pathLooksUsable(socketPath),
  };
}

// Resolve Podman API endpoint: info → machine pipe/socket → system connection.
function resolvePodmanOptionsFromCli(podmanBin) {
  const info = readPodmanInfo(podmanBin);
  const remote = remoteSocketFromInfo(info);
  const rootless = isRootlessFromPodmanInfo(info) || !isRootUser();

  if (remote.path && remote.exists) {
    return {
      options: { socketPath: remote.path },
      rootless,
    };
  }

  // Windows/macOS: Podman runs inside a VM; the host talks over a pipe/socket.
  try {
    const machines = parseJsonSafe(runCliSync(podmanBin, ["machine", "inspect"]));
    const machine = Array.isArray(machines) ? machines[0] : machines;
    const pipePath = machine?.ConnectionInfo?.PodmanPipe?.Path;
    const socketPath = machine?.ConnectionInfo?.PodmanSocket?.Path;
    if (process.platform === "win32" && pipePath) {
      return { options: parseHostToOptions(pipePath), rootless: false, info };
    }
    if (socketPath && pathLooksUsable(socketPath.replace(/^unix:\/\//, ""))) {
      return { options: parseHostToOptions(socketPath), rootless: false, info };
    }
  } catch {
    // machine inspect is only available when a VM-backed install is present.
  }

  try {
    const connections = parseJsonSafe(
      runCliSync(podmanBin, ["system", "connection", "list", "--format", "json"])
    );
    const list = Array.isArray(connections) ? connections : [];
    const selected = list.find((item) => item.Default) || list[0];
    if (selected?.URI) {
      const options = parseHostToOptions(selected.URI);
      if (options.host || pathLooksUsable(options.socketPath)) {
        return { options, rootless, info };
      }
    }
  } catch {
    // Older Podman builds may not support system connection list.
  }

  return {
    options: null,
    rootless,
    intendedSocket: remote.path || preferredLocalPodmanSocket(rootless),
    socketMissing: true,
    info,
  };
}

// Normalize a found engine into one object the rest of the app can store/use.
function describeRuntime(engine, options, extra = {}) {
  const socketPath = options?.socketPath || "";
  const rootless =
    extra.rootless === true ||
    (engine === "podman" &&
      process.platform === "linux" &&
      (socketPath.includes(`${path.sep}run${path.sep}user${path.sep}`) ||
        socketPath.includes("/run/user/")));

  return {
    engine,
    options,
    socketPath,
    rootless,
    binaryPath: extra.binaryPath || null,
    detected: Boolean(options && (options.socketPath || options.host)),
  };
}

// Locate Docker without starting anything. Returns null if not found.
function tryResolveDocker() {
  const dockerBin = findExecutable("docker");
  const envOptions = optionsFromEnv("docker");

  let options = envOptions;
  if (!options && dockerBin) {
    try {
      options = resolveDockerOptionsFromCli(dockerBin);
    } catch (error) {
      if (error?.code !== "ENOENT") {
        // Context inspect can fail while Docker Desktop is still booting.
        options = null;
      }
    }
  }

  const preferred = options?.socketPath;
  const socketPath = firstUsablePath(dockerSocketCandidates(preferred), preferred);
  if (socketPath) {
    return describeRuntime("docker", { ...(options || {}), socketPath }, { binaryPath: dockerBin });
  }
  if (options?.host) {
    return describeRuntime("docker", options, { binaryPath: dockerBin });
  }
  if (dockerBin && process.platform === "win32") {
    return describeRuntime(
      "docker",
      { socketPath: "\\\\.\\pipe\\dockerDesktopLinuxEngine" },
      { binaryPath: dockerBin }
    );
  }
  if (dockerBin) {
    return describeRuntime("docker", { socketPath: "/var/run/docker.sock" }, { binaryPath: dockerBin });
  }
  return null;
}

// Locate Podman without starting the machine/socket. Returns null if not found.
function tryResolvePodman() {
  const podmanBin = findExecutable("podman");
  const envOptions = optionsFromEnv("podman");
  if (envOptions?.host) {
    return withPodmanMachineMeta(
      describeRuntime("podman", envOptions, { binaryPath: podmanBin, rootless: false }),
      podmanBin
    );
  }
  if (envOptions?.socketPath && pathLooksUsable(envOptions.socketPath)) {
    return withPodmanMachineMeta(
      describeRuntime("podman", envOptions, { binaryPath: podmanBin, rootless: !isRootUser() }),
      podmanBin
    );
  }

  let cliResult = { options: null, rootless: !isRootUser() };
  if (podmanBin) {
    try {
      cliResult = resolvePodmanOptionsFromCli(podmanBin);
    } catch {
      cliResult = { options: null, rootless: !isRootUser() };
    }
  }

  const preferred = cliResult.options?.socketPath;
  const socketPath = firstUsablePath(podmanSocketCandidates(preferred), null);
  if (socketPath) {
    return withPodmanMachineMeta(
      describeRuntime(
        "podman",
        { ...(cliResult.options || {}), socketPath },
        { binaryPath: podmanBin, rootless: cliResult.rootless }
      ),
      podmanBin
    );
  }
  if (cliResult.options?.host) {
    return withPodmanMachineMeta(
      describeRuntime("podman", cliResult.options, {
        binaryPath: podmanBin,
        rootless: cliResult.rootless,
      }),
      podmanBin
    );
  }
  if (
    cliResult.options?.socketPath &&
    process.platform === "win32" &&
    cliResult.options.socketPath.includes("\\pipe\\")
  ) {
    return withPodmanMachineMeta(
      describeRuntime("podman", cliResult.options, {
        binaryPath: podmanBin,
        rootless: false,
      }),
      podmanBin
    );
  }
  return null;
}

function findPodmanBinary() {
  return findExecutable("podman");
}

// Linux: start the podman.socket systemd unit so the API socket appears.
async function startPodmanSocketUnit(rootless) {
  const systemctl = findExecutable("systemctl");
  if (!systemctl) {
    return;
  }
  try {
    if (rootless && !isRootUser()) {
      await runCli(systemctl, ["--user", "start", "podman.socket"], 5000);
      return;
    }
    await runCli(systemctl, ["start", "podman.socket"], 5000);
  } catch {
    // Unit may not be installed (macOS, or distros without the socket unit).
  }
}

// Last resort: run Podman's own HTTP API server on the expected socket path.
function startPodmanSystemService(podmanBin, socketPath) {
  if (!podmanBin || !socketPath) {
    return;
  }
  if (podmanServiceChild && !podmanServiceChild.killed) {
    return;
  }
  try {
    fs.mkdirSync(path.dirname(socketPath), { recursive: true });
  } catch {
    // Directory may already exist or be created by podman itself.
  }
  const uri = socketPath.startsWith("unix://") ? socketPath : `unix://${socketPath}`;
  podmanServiceChild = spawn(podmanBin, ["system", "service", "--time=0", uri], {
    detached: true,
    stdio: "ignore",
    windowsHide: true,
  });
  podmanServiceChild.unref();
}

// Windows/macOS: start the Podman VM if it exists but is stopped.
async function ensurePodmanMachineRunning(podmanBin) {
  try {
    const listRaw = await runCli(podmanBin, ["machine", "list", "--format", "json"], CLI_TIMEOUT_MS);
    const machines = parseMachineList(listRaw);
    if (!machines.length) {
      return false;
    }
    const machineName = defaultMachineName(machines);
    const machine =
      machines.find((item) => item.Name === machineName) || machines[0];
    const running =
      machine.Running === true ||
      machine.Running === "true" ||
      String(machine.LastUp || "").toLowerCase().includes("currently running");
    if (!running) {
      const nameArgs = machineName ? [machineName] : [];
      await runCli(podmanBin, ["machine", "start", ...nameArgs], MACHINE_TIMEOUT_MS);
    }

    // After start (or WSL reboot), the named pipe/socket can lag a few seconds
    // behind "machine started". Wait so probes don't hit connect ENOENT.
    if (process.platform === "win32") {
      const pipePath = `\\\\.\\pipe\\${machineName || "podman-machine-default"}`;
      await waitForPath(pipePath, 20000);
    } else {
      const after = tryResolvePodman();
      if (after?.socketPath) {
        await waitForPath(after.socketPath, 20000);
      }
    }
    return true;
  } catch (error) {
    console.warn("ensurePodmanMachineRunning failed:", error?.message || error);
    return false;
  }
}

// Make Podman usable: start machine (Win/mac) or socket/service (Linux) if needed.
async function ensurePodmanRuntime() {
  const podmanBin = findPodmanBinary();
  if (!podmanBin) {
    return null;
  }

  const already = tryResolvePodman();
  if (already) {
    return already;
  }

  if (process.platform === "darwin" || process.platform === "win32") {
    await ensurePodmanMachineRunning(podmanBin);
    const afterMachine = tryResolvePodman();
    if (afterMachine) {
      return afterMachine;
    }
  }

  const info = readPodmanInfo(podmanBin);
  const remote = remoteSocketFromInfo(info);
  const rootless = isRootlessFromPodmanInfo(info) || !isRootUser();
  const rootfulSockets = new Set([defaultRootfulPodmanSocket(), "/var/run/podman/podman.sock"]);
  let socketPath = remote.path;
  if (rootless && rootfulSockets.has(socketPath)) {
    socketPath = preferredLocalPodmanSocket(true);
  }
  socketPath = socketPath || preferredLocalPodmanSocket(rootless);

  if (pathLooksUsable(socketPath)) {
    return withPodmanMachineMeta(
      describeRuntime("podman", { socketPath }, { binaryPath: podmanBin, rootless }),
      podmanBin
    );
  }

  await startPodmanSocketUnit(rootless);
  if (await waitForPath(socketPath, 4000)) {
    return withPodmanMachineMeta(
      describeRuntime("podman", { socketPath }, { binaryPath: podmanBin, rootless }),
      podmanBin
    );
  }

  startPodmanSystemService(podmanBin, socketPath);
  if (await waitForPath(socketPath, 8000)) {
    return withPodmanMachineMeta(
      describeRuntime("podman", { socketPath }, { binaryPath: podmanBin, rootless }),
      podmanBin
    );
  }

  return tryResolvePodman();
}

// Lightweight "is Docker / Podman installed?" snapshot for the UI.
function inspectAvailableRuntimes() {
  const docker = tryResolveDocker();
  const podmanBin = findPodmanBinary();
  const podman = tryResolvePodman();
  return {
    docker: {
      detected: Boolean(docker || findExecutable("docker")),
      socketPath: docker?.socketPath || "",
      binaryPath: docker?.binaryPath || findExecutable("docker"),
    },
    podman: {
      detected: Boolean(podmanBin),
      socketPath: podman?.socketPath || "",
      binaryPath: podmanBin,
      rootless: Boolean(podman?.rootless || (podmanBin && !isRootUser())),
    },
  };
}

// Engines to try, in order. Auto prefers Docker when both are available.
function listRuntimeCandidates(preference = readRuntimePreference()) {
  const docker = tryResolveDocker();
  const podman = tryResolvePodman();

  if (preference === "docker") {
    return docker ? [docker] : [];
  }
  if (preference === "podman") {
    return podman ? [podman] : [];
  }
  return [docker, podman].filter(Boolean);
}

// Sync pick of the preferred/available engine (used by $docker getter).
function resolveContainerRuntimeSync(forceRefresh = false) {
  if (!forceRefresh && cachedRuntime) {
    return cachedRuntime;
  }

  const preference = readRuntimePreference();
  const candidates = listRuntimeCandidates(preference);
  if (candidates.length === 0) {
    const wanted = preference === "auto" ? "Docker or Podman" : preference === "podman" ? "Podman" : "Docker";
    throw new Error(`${wanted} is not installed or its socket could not be located.`);
  }

  cachedRuntime = candidates[0];
  return cachedRuntime;
}

function setCachedRuntime(runtime) {
  cachedRuntime = runtime || null;
}

function clearRuntimeCache() {
  cachedRuntime = null;
}

function setRuntimePreference(preference) {
  persistRuntimePreference(preference);
  clearRuntimeCache();
  return preference;
}

function getCachedRuntime() {
  return cachedRuntime;
}

function getDockerodeOptionsFromContextSync() {
  try {
    return resolveContainerRuntimeSync().options;
  } catch (error) {
    if (error?.code === "ENOENT") {
      const fallbackSocketPath = firstUsablePath(dockerSocketCandidates(), null);
      if (fallbackSocketPath) {
        return { socketPath: fallbackSocketPath };
      }
      throw new Error("Docker CLI is not installed or not available in PATH.");
    }
    throw error;
  }
}

// Options object for `new Docker(...)` — always re-reads cache so preference switches apply.
function getResolvedDockerodeOptions() {
  if (cachedRuntime?.options) {
    return cachedRuntime.options;
  }
  try {
    return resolveContainerRuntimeSync().options;
  } catch {
    return {};
  }
}

// Podman's /version response usually contains "podman"; Docker's does not.
function identifyEngineFromVersion(versionInfo) {
  const blob = JSON.stringify(versionInfo || {}).toLowerCase();
  if (blob.includes("podman")) {
    return "podman";
  }
  return "docker";
}

// True when the last ":..." part of a bind is options (ro/rw/Z), not a path.
function lastSegmentIsMountOptions(segment) {
  if (!segment) {
    return false;
  }
  if (segment.startsWith("/") || segment.startsWith("\\")) {
    return false;
  }
  if (/^[A-Za-z]$/.test(segment)) {
    return false;
  }
  return /^(ro|rw|z|Z)(,|$)/.test(segment) || segment.includes(",");
}

// Split "host:container[:opts]" safely — Windows drives use "C:" too.
function parseBindSpec(bind) {
  if (!bind || typeof bind !== "string") {
    return null;
  }
  const windowsHost = bind.match(/^([A-Za-z]:[\\/][^:]*):(.*)$/);
  let host;
  let remainder;
  if (windowsHost) {
    host = windowsHost[1];
    remainder = windowsHost[2];
  } else {
    const idx = bind.indexOf(":");
    if (idx === -1) {
      return null;
    }
    host = bind.slice(0, idx);
    remainder = bind.slice(idx + 1);
  }

  const optionSep = remainder.lastIndexOf(":");
  if (optionSep > 0 && lastSegmentIsMountOptions(remainder.slice(optionSep + 1))) {
    return {
      host,
      container: remainder.slice(0, optionSep),
      options: remainder.slice(optionSep + 1),
    };
  }
  return { host, container: remainder, options: "" };
}

function formatBindSpec({ host, container, options }) {
  if (!host || !container) {
    return null;
  }
  return options ? `${host}:${container}:${options}` : `${host}:${container}`;
}

// C:\Users\... → /mnt/c/Users/... (how Podman's Windows VM sees the host disk).
function windowsPathToWslMount(hostPath) {
  const normalized = String(hostPath || "").replace(/\\/g, "/");
  const match = normalized.match(/^([A-Za-z]):\/(.*)$/);
  if (!match) {
    return normalized;
  }
  return `/mnt/${match[1].toLowerCase()}/${match[2]}`;
}

// Rewrite Windows host paths for Podman: use machine mounts if known, else /mnt/<drive>/.
function rewriteHostPathForPodman(hostPath, runtime) {
  if (process.platform !== "win32") {
    return String(hostPath || "").replace(/\\/g, "/");
  }
  const normalized = String(hostPath || "").replace(/\\/g, "/");
  const mounts = Array.isArray(runtime?.machineMounts) ? runtime.machineMounts : [];
  for (const mount of mounts) {
    const source = String(mount.Source || mount.source || "").replace(/\\/g, "/");
    const destination = mount.Destination || mount.destination || mount.Target || mount.target;
    if (!source || !destination) {
      continue;
    }
    if (normalized.toLowerCase().startsWith(source.toLowerCase())) {
      const suffix = normalized.slice(source.length);
      return `${destination}${suffix}`.replace(/\\/g, "/");
    }
  }
  return windowsPathToWslMount(normalized);
}

// Read VM type + shared folders from `podman machine inspect` (Win/mac only).
function inspectPodmanMachine(podmanBin) {
  if (!podmanBin || (process.platform !== "win32" && process.platform !== "darwin")) {
    return { vmType: "", mounts: [] };
  }
  try {
    const machines = parseJsonSafe(runCliSync(podmanBin, ["machine", "inspect"]));
    const machine = Array.isArray(machines) ? machines[0] : machines;
    if (!machine) {
      return { vmType: "", mounts: [] };
    }
    return {
      vmType: String(machine.VMType || machine.Provider || "").toLowerCase(),
      mounts: Array.isArray(machine.Mounts) ? machine.Mounts : [],
    };
  } catch {
    return { vmType: "", mounts: [] };
  }
}

function withPodmanMachineMeta(runtime, podmanBin) {
  if (!runtime) {
    return runtime;
  }
  const meta = inspectPodmanMachine(podmanBin);
  return {
    ...runtime,
    vmType: meta.vmType || runtime.vmType || "",
    machineMounts: meta.mounts.length ? meta.mounts : runtime.machineMounts || [],
  };
}

// Append a mount option like :Z without duplicating it.
function appendBindOption(bind, option) {
  if (!bind) {
    return bind;
  }
  const optionPattern = new RegExp(`(^|,)${option}(,|$)`);
  const lastColon = bind.lastIndexOf(":");
  if (lastColon === -1) {
    return `${bind}:${option}`;
  }
  const lastPart = bind.slice(lastColon + 1);
  if (lastSegmentIsMountOptions(lastPart)) {
    if (optionPattern.test(lastPart) || /(^|,)z(,|$)/i.test(lastPart)) {
      return bind;
    }
    return `${bind},${option}`;
  }
  return `${bind}:${option}`;
}

/**
 * Final pass on HostConfig.Binds before create/run:
 * - Windows Podman: rewrite C:\... → /mnt/c/...
 * - Linux Podman: add :Z for SELinux
 */
function prepareBindMounts(binds) {
  if (!Array.isArray(binds)) {
    return binds;
  }
  const runtime = cachedRuntime || {};
  return binds.map((bind) => {
    const parsed = parseBindSpec(bind);
    if (!parsed) {
      return bind;
    }
    let { host, container, options } = parsed;
    if (runtime.engine === "podman" && process.platform === "win32") {
      host = rewriteHostPathForPodman(host, runtime);
    }
    let formatted = formatBindSpec({ host, container, options }) || bind;
    if (runtime.engine === "podman" && process.platform === "linux") {
      formatted = appendBindOption(formatted, "Z");
    }
    return formatted;
  });
}

/**
 * Which user the container process should run as.
 * Podman rootless / Win / mac usually need 0:0 so bind mounts are writable.
 */
function getContainerUser(hostUid, hostGid) {
  const runtime = cachedRuntime || null;
  if (runtime?.engine === "podman") {
    if (process.platform === "linux" && runtime.rootless) {
      return "0:0";
    }
    if (process.platform === "win32" || process.platform === "darwin") {
      return "0:0";
    }
  }
  if (hostUid == null || hostGid == null) {
    return undefined;
  }
  return `${hostUid}:${hostGid}`;
}

function engineDisplayName(engine) {
  if (engine === "podman") {
    return "Podman";
  }
  if (engine === "docker") {
    return "Docker";
  }
  return "container engine";
}

function parseMachineList(raw) {
  const parsed = parseJsonSafe(raw);
  if (Array.isArray(parsed)) {
    return parsed;
  }
  if (parsed && Array.isArray(parsed.list)) {
    return parsed.list;
  }
  return [];
}

function defaultMachineName(machines) {
  if (!machines.length) {
    return "";
  }
  const selected =
    machines.find((machine) => machine.Default || machine.Default === "true") ||
    machines.find((machine) => machine.Running || machine.Running === "true") ||
    machines[0];
  return selected?.Name || "";
}

async function getPodmanBinary() {
  const cached = cachedRuntime?.engine === "podman" ? cachedRuntime.binaryPath : null;
  return cached || findExecutable("podman");
}

// Win/mac Resource Manager: apply CPU/RAM then restart the machine.
// On Windows, Podman uses WSL — `podman machine set --cpus/--memory` is not
// supported there. Resources are global for all WSL2 distros via ~/.wslconfig.
async function applyPodmanMachineResources({ cpus, memoryMiB }) {
  const podmanBin = await getPodmanBinary();
  if (!podmanBin) {
    throw new Error("Podman CLI is not installed or not available in PATH.");
  }

  const listRaw = await runCli(podmanBin, ["machine", "list", "--format", "json"], CLI_TIMEOUT_MS);
  const machineName = defaultMachineName(parseMachineList(listRaw));
  const nameArgs = machineName ? [machineName] : [];

  try {
    await runCli(podmanBin, ["machine", "stop", ...nameArgs], MACHINE_TIMEOUT_MS);
  } catch (error) {
    const message = String(error?.stderr || error?.message || "");
    if (!/not running|already stopped/i.test(message)) {
      throw error;
    }
  }

  if (process.platform === "win32") {
    writeWslResourceConfig({
      memoryGb: Math.max(1, Math.round(Number(memoryMiB) / 1024)),
      processors: Math.max(1, Number(cpus) || 1),
    });
    // WSL must fully shut down before .wslconfig is picked up.
    try {
      await runCli("wsl.exe", ["--shutdown"], CLI_TIMEOUT_MS);
    } catch (error) {
      const message = String(error?.stderr || error?.message || "");
      if (!/no running distributions|there are no distributions/i.test(message)) {
        throw error;
      }
    }
    await new Promise((resolve) => setTimeout(resolve, 8000));
  } else {
    // QEMU-backed machines (typical on some Linux/mac setups) support per-VM limits.
    await runCli(
      podmanBin,
      ["machine", "set", "--cpus", String(cpus), "--memory", String(memoryMiB), ...nameArgs],
      CLI_TIMEOUT_MS
    );
  }

  await runCli(podmanBin, ["machine", "start", ...nameArgs], MACHINE_TIMEOUT_MS);
  clearRuntimeCache();
  return {
    machineName: machineName || "default",
    appliedVia: process.platform === "win32" ? "wslconfig" : "podman-machine-set",
  };
}

function writeWslResourceConfig({ memoryGb, processors }) {
  const wslConfigPath = path.join(os.homedir(), ".wslconfig");
  let content = "";
  if (fs.existsSync(wslConfigPath)) {
    content = fs.readFileSync(wslConfigPath, "utf8");
  }
  if (!/\[wsl2\]/i.test(content)) {
    content = content.trim().length ? `${content.trim()}\n\n[wsl2]\n` : "[wsl2]\n";
  }
  if (/^\s*memory\s*=/im.test(content)) {
    content = content.replace(/^\s*memory\s*=\s*[^\r\n]+/im, `memory=${memoryGb}GB`);
  } else {
    content = content.replace(/\[wsl2\]/i, `[wsl2]\nmemory=${memoryGb}GB`);
  }
  if (/^\s*processors\s*=/im.test(content)) {
    content = content.replace(/^\s*processors\s*=\s*[^\r\n]+/im, `processors=${processors}`);
  } else {
    content = content.replace(/\[wsl2\]/i, `[wsl2]\nprocessors=${processors}`);
  }
  fs.writeFileSync(wslConfigPath, content, "utf8");
}

module.exports = {
  PREFERENCES,
  appendBindOption,
  applyPodmanMachineResources,
  clearRuntimeCache,
  engineDisplayName,
  findExecutable,
  getCachedRuntime,
  getContainerUser,
  getDockerodeOptionsFromContextSync,
  getResolvedDockerodeOptions,
  identifyEngineFromVersion,
  inspectAvailableRuntimes,
  listRuntimeCandidates,
  ensurePodmanRuntime,
  parseHostToOptions,
  prepareBindMounts,
  readRuntimePreference,
  resolveContainerRuntimeSync,
  setCachedRuntime,
  setRuntimePreference,
};
