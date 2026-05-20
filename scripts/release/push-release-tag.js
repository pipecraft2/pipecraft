/**
 * Push an annotated-style release tag for a remote branch using version from that branch's package.json.
 * Triggers .github/workflows/release.yml (tag pattern v*).
 *
 * Usage:
 *   node scripts/release/push-release-tag.js <branch> [--remote fork]
 *   yarn release:push -- main
 *   yarn release:push -- test-v1.2.0 --remote fork
 */
/* eslint-disable no-console */
const { spawnSync } = require("child_process");

function git(args) {
  return spawnSync("git", args, { encoding: "utf8", shell: false });
}

function die(msg) {
  console.error(msg);
  process.exit(1);
}

function parseArgs(argv) {
  let remote = null;
  let branch = null;
  for (let i = 0; i < argv.length; i += 1) {
    const a = argv[i];
    if (a === "--remote") {
      remote = argv[++i];
      if (!remote) die("--remote needs a value");
    } else if (!a.startsWith("-") && branch === null) {
      branch = a;
    } else if (a.startsWith("-")) {
      die(`Unknown option: ${a}`);
    }
  }
  return { remote, branch };
}

function defaultRemote() {
  const r = git(["remote", "get-url", "origin"]);
  return r.status === 0 ? "origin" : "fork";
}

function main() {
  const { remote: remoteArg, branch } = parseArgs(process.argv.slice(2));
  if (!branch) {
    die(
      "Usage: node scripts/release/push-release-tag.js <branch> [--remote fork]\n" +
        "Example: yarn release:push -- main"
    );
  }

  const remote = remoteArg || defaultRemote();

  let r = git(["fetch", remote, branch]);
  if (r.status !== 0) {
    die(`git fetch ${remote} ${branch} failed:\n${r.stderr || r.stdout}`);
  }

  const ref = `${remote}/${branch}`;
  r = git(["rev-parse", ref]);
  if (r.status !== 0) {
    die(`git rev-parse ${ref} failed:\n${r.stderr || r.stdout}`);
  }
  const sha = r.stdout.trim();

  r = git(["show", `${sha}:package.json`]);
  if (r.status !== 0) {
    die(`Cannot read package.json at ${sha}:\n${r.stderr || r.stdout}`);
  }

  let version;
  try {
    version = JSON.parse(r.stdout).version;
  } catch (e) {
    die(`Invalid package.json at ${sha}: ${e.message}`);
  }
  if (!version || typeof version !== "string") {
    die(`Missing or invalid "version" in package.json at ${ref} (${sha})`);
  }

  const tag = `v${version}`;

  r = git(["ls-remote", remote, `refs/tags/${tag}`]);
  if (r.status !== 0) {
    die(`git ls-remote failed:\n${r.stderr || r.stdout}`);
  }
  if (r.stdout.trim()) {
    die(
      `Tag ${tag} already exists on ${remote}. Bump "version" in package.json on ${branch}, or delete the remote tag.`
    );
  }

  r = git(["push", remote, `${sha}:refs/tags/${tag}`]);
  if (r.status !== 0) {
    die(`git push failed:\n${r.stderr || r.stdout}`);
  }

  console.log(`Pushed tag ${tag} (${sha}) from ${remote}/${branch}.`);
  console.log('GitHub Actions workflow "Release" should start for this tag.');
}

main();
