<template>
  <v-container fluid class="rm-page">
    <v-card dark class="rm-card">
      <div class="rm-header">
        <span class="rm-title">RESOURCE MANAGER</span>
        <v-spacer></v-spacer>
        <v-icon x-small :color="statusColor" class="mr-2">mdi-circle</v-icon>
        <span class="rm-status" role="status">{{ runtimeStatusLine }}</span>
      </div>
      <v-divider></v-divider>

      <v-row no-gutters class="rm-columns">
        <v-col cols="12" sm="6" :class="['rm-section', { 'rm-col-divider': $vuetify.breakpoint.smAndUp }]">
          <div class="rm-section-title">Engine</div>
          <template v-if="bothEnginesDetected">
            <v-btn-toggle
              :value="selectedEngine"
              :mandatory="Boolean(selectedEngine)"
              dense
              color="teal accent-3"
              background-color="transparent"
              :disabled="runtimeStarting"
              @change="handleEngineChange"
            >
              <v-btn value="docker" small :disabled="!availableRuntimes.docker">Docker</v-btn>
              <v-btn value="podman" small :disabled="!availableRuntimes.podman">Podman</v-btn>
            </v-btn-toggle>
            <div class="d-flex align-center mt-2">
              <v-checkbox
                dense
                hide-details
                class="mt-0 pt-0"
                color="teal accent-3"
                label="Use as default"
                :disabled="!selectedEngine || runtimeStarting"
                :input-value="useAsDefault"
                @change="handleDefaultChange"
              ></v-checkbox>
              <v-tooltip bottom max-width="260">
                <template v-slot:activator="{ on, attrs }">
                  <v-icon small class="ml-2 rm-muted" aria-label="About Use as default" v-bind="attrs" v-on="on">
                    mdi-information-outline
                  </v-icon>
                </template>
                When checked, this engine starts automatically next time. Uncheck to choose again at launch.
              </v-tooltip>
            </div>
          </template>
          <div class="rm-muted mt-2">Detected: {{ detectedRuntimesText }}</div>
          <div v-if="socketPath" class="rm-muted rm-break">Socket: {{ socketPath }}</div>
          <v-btn
            v-if="startAction"
            small
            depressed
            color="teal accent-3"
            class="black--text mt-3"
            :loading="runtimeStarting"
            :disabled="runtimeStarting"
            @click="startAction.handler"
          >
            {{ startAction.label }}
          </v-btn>
        </v-col>

        <v-col cols="12" sm="6" class="rm-section">
          <div class="rm-section-title">Resources per workflow</div>
          <div class="rm-slider-row">
            <span class="rm-slider-label">CPU</span>
            <v-slider
              v-model="ncpu"
              dense
              hide-details
              color="teal accent-3"
              track-color="grey darken-1"
              step="1"
              min="1"
              :max="processors"
              :disabled="!isDockerActive"
              aria-label="CPU cores per workflow"
            ></v-slider>
            <span class="rm-slider-value">{{ isDockerActive ? `${ncpu} / ${processors}` : "–" }}</span>
          </div>
          <div class="rm-slider-row">
            <span class="rm-slider-label">RAM</span>
            <v-slider
              v-model="memtotal"
              dense
              hide-details
              color="teal accent-3"
              track-color="grey darken-1"
              step="1"
              min="1"
              :max="memory"
              :disabled="!isDockerActive"
              aria-label="Memory in GB per workflow"
            ></v-slider>
            <span class="rm-slider-value">{{ isDockerActive ? `${memtotal} / ${memory} GB` : "–" }}</span>
          </div>
          <div class="d-flex align-center mt-2">
            <span class="rm-muted">{{ resourceHint }}</span>
            <v-spacer></v-spacer>
            <v-btn
              v-if="isDockerActive && showApplyRestart"
              small
              depressed
              color="teal accent-3"
              class="black--text ml-2"
              :loading="runtimeStarting"
              :disabled="runtimeStarting"
              @click="handleApplyResources"
            >
              Apply &amp; restart
            </v-btn>
          </div>
        </v-col>
      </v-row>
      <v-divider></v-divider>

      <div class="rm-section">
        <div class="d-flex align-center">
          <span class="rm-section-title mb-0">Database cache</span>
          <v-chip x-small class="ml-2">{{ cacheEntries.length }} · {{ formatBytes(cacheTotalBytes) }}</v-chip>
          <v-tooltip bottom max-width="300">
            <template v-slot:activator="{ on, attrs }">
              <v-icon small class="ml-2 rm-muted" aria-label="About the database cache" v-bind="attrs" v-on="on">
                mdi-help-circle-outline
              </v-icon>
            </template>
            Search indexes built from reference databases (BLAST, SINTAX, NextITS chimera database),
            reused by every project. Removing an entry is safe; it is rebuilt the next time a workflow needs it.
          </v-tooltip>
          <v-spacer></v-spacer>
          <v-btn icon small title="Refresh" aria-label="Refresh database cache list" @click="handleRefreshCache">
            <v-icon small>mdi-refresh</v-icon>
          </v-btn>
          <v-btn icon small title="Open folder" aria-label="Open database cache folder" @click="handleOpenCacheFolder">
            <v-icon small>mdi-folder-open-outline</v-icon>
          </v-btn>
        </div>
        <v-data-table
          v-model="selectedCacheEntries"
          :headers="cacheHeaders"
          :items="cacheEntries"
          item-key="id"
          show-select
          checkbox-color="teal accent-3"
          dense
          disable-pagination
          hide-default-footer
          no-data-text="The database cache is empty."
          class="rm-table mt-1"
        >
          <template v-slot:[`item.name`]="{ item }">
            <span :title="item.sourcePath">{{ item.name }}</span>
          </template>
          <template v-slot:[`item.sizeBytes`]="{ item }">
            {{ formatBytes(item.sizeBytes) }}
          </template>
          <template v-slot:[`item.lastUsedAt`]="{ item }">
            <span :title="formatDateTime(item.lastUsedAt)">{{ formatRelativeDate(item.lastUsedAt) }}</span>
          </template>
          <template v-slot:[`item.status`]="{ item }">
            <v-chip
              x-small
              label
              outlined
              :color="cacheStatusInfo(item.status).color"
              :title="cacheStatusInfo(item.status).hint"
            >
              {{ cacheStatusInfo(item.status).label }}
            </v-chip>
          </template>
        </v-data-table>
        <div class="d-flex flex-wrap align-center justify-end mt-2">
          <v-btn
            small
            outlined
            class="ml-2"
            :disabled="cacheLocked || !selectedCacheEntries.length"
            @click="handleRemoveSelected"
          >
            Remove selected{{ selectedCacheEntries.length ? ` (${selectedCacheEntries.length})` : "" }}
          </v-btn>
          <v-btn
            small
            outlined
            class="ml-2"
            :disabled="cacheLocked || !outdatedCacheEntries.length"
            @click="handleRemoveOutdated"
          >
            Remove outdated{{ outdatedCacheEntries.length ? ` (${outdatedCacheEntries.length})` : "" }}
          </v-btn>
          <v-btn
            small
            outlined
            color="red lighten-2"
            class="ml-2"
            :disabled="cacheLocked || !cacheEntries.length"
            @click="handleRemoveAll"
          >
            Remove all
          </v-btn>
        </div>
        <div v-if="cacheLocked" class="rm-muted text-right mt-1">
          Cache entries can be removed after the running workflow finishes.
        </div>
      </div>
    </v-card>
  </v-container>
</template>

<script>
const { exec, execFile } = require("child_process");
const fs = require("fs");
import os from "os";
import { mapState, mapGetters } from "vuex";
import { applyPodmanMachineResources, dockerDesktopWindowsExe, writeWslResourceConfig } from "../utils/containerRuntime";
import { getCacheRoot, listCacheEntries, removeCacheEntries } from "../utils/dbCache";
const CPU = os.cpus().length;
const MEM = Number((os.totalmem() / 1024 ** 3).toFixed(0));
const Swal = require("sweetalert2");

const CACHE_STATUS = {
  ok: { label: "Ready", color: "teal accent-3", hint: "Up to date with the source file." },
  incomplete: { label: "Incomplete", color: "grey lighten-1", hint: "The build did not finish; it is rebuilt when needed." },
  missing: { label: "Source missing", color: "orange lighten-2", hint: "The source file was moved or deleted." },
  changed: { label: "Source changed", color: "orange lighten-2", hint: "The source file was edited; a new index is built for the new version." },
};
const OUTDATED_STATUSES = ["incomplete", "missing", "changed"];

const DAY_MS = 24 * 60 * 60 * 1000;

export default {
  name: "ResourceManager",
  computed: {
    ...mapState({
      dockerSettingsPath: state => state.systemSpecs.dockerSettings,
      runtimePreference: state => state.containerRuntime.preference,
      runtimeStarting: state => state.containerRuntime.starting,
      socketPath: state => state.containerRuntime.socketPath,
      availableRuntimes: state => state.containerRuntime.available,
    }),
    ...mapGetters(['isDockerActive', 'engineLabel', 'runtimeStatusText', 'activeEngine', 'bothEnginesDetected', 'sessionEngine']),
    selectedEngine() {
      return this.activeEngine || this.sessionEngine || this.runtimePreference || "";
    },
    useAsDefault() {
      return Boolean(this.selectedEngine && this.runtimePreference === this.selectedEngine);
    },
    ncpu: {
      get() {
        return this.$store.state.dockerInfo.NCPU;
      },
      set(value) {
        this.$store.commit("setNCPU", value);
      },
    },
    memtotal: {
      get() {
        return Math.ceil(
          (this.$store.state.dockerInfo.MemTotal / 1024 ** 3).toFixed(2)
        );
      },
      set(value) {
        this.$store.commit("setMemTotal", value * 1024 ** 3);
      },
    },
    resourceEngine() {
      return this.selectedEngine;
    },
    showApplyRestart() {
      return this.$store.state.OStype !== "Linux";
    },
    resourceHint() {
      if (!this.isDockerActive) return "Start an engine to change limits.";
      if (!this.showApplyRestart) return "Applied to each workflow container. No restart needed.";
      return this.resourceEngine === "podman"
        ? "Restarts the Podman machine to apply."
        : "Restarts Docker Desktop to apply.";
    },
    startAction() {
      if (this.isDockerActive) return null;
      if (!this.selectedEngine) {
        return { label: "Check for Docker or Podman", handler: this.handleDetectEngines };
      }
      return { label: `Start ${this.engineLabel}`, handler: this.handleStartEngine };
    },
    statusColor() {
      if (this.runtimeStarting) return "amber";
      return this.isDockerActive ? "teal accent-3" : "grey";
    },
    runtimeStatusLine() {
      const rootless = this.$store.state.containerRuntime.rootless ? " (rootless)" : "";
      return `${this.runtimeStatusText}${this.isDockerActive ? rootless : ""}`;
    },
    detectedRuntimesText() {
      const found = [];
      if (this.availableRuntimes.docker) found.push("Docker");
      if (this.availableRuntimes.podman) found.push("Podman");
      return found.length ? found.join(", ") : "none";
    },
    // A running workflow may be building or reading a cache entry.
    cacheLocked() {
      return Boolean(this.$store.state.runInfo.active);
    },
    cacheTotalBytes() {
      return this.cacheEntries.reduce((total, entry) => total + entry.sizeBytes, 0);
    },
    outdatedCacheEntries() {
      return this.cacheEntries.filter((entry) => OUTDATED_STATUSES.includes(entry.status));
    },
  },
  data() {
    return {
      processors: CPU,
      memory: MEM,
      cacheEntries: [],
      selectedCacheEntries: [],
      cacheHeaders: [
        { text: "Database", value: "name" },
        { text: "Type", value: "format" },
        { text: "Size", value: "sizeBytes" },
        { text: "Last used", value: "lastUsedAt" },
        { text: "Status", value: "status" },
      ],
    };
  },
  watch: {
    cacheLocked(locked) {
      if (!locked) this.handleRefreshCache();
    },
  },
  mounted() {
    this.handleRefreshCache();
  },
  methods: {
    handleRefreshCache() {
      try {
        this.cacheEntries = listCacheEntries();
      } catch (error) {
        console.error("Reading the database cache failed:", error);
        this.cacheEntries = [];
      }
      const ids = new Set(this.cacheEntries.map((entry) => entry.id));
      this.selectedCacheEntries = this.selectedCacheEntries.filter((entry) => ids.has(entry.id));
    },
    async handleOpenCacheFolder() {
      const root = getCacheRoot();
      try {
        fs.mkdirSync(root, { recursive: true });
        const error = await require("@electron/remote").shell.openPath(root);
        if (error) throw new Error(error);
      } catch (error) {
        console.error(`Opening ${root} failed:`, error);
      }
    },
    handleRemoveSelected() {
      this.removeCache(this.selectedCacheEntries);
    },
    handleRemoveOutdated() {
      this.removeCache(this.outdatedCacheEntries);
    },
    async handleRemoveAll() {
      const confirmation = await Swal.fire({
        title: "Remove all cached databases?",
        text: `This frees ${this.formatBytes(this.cacheTotalBytes)}. Indexes are rebuilt the next time a workflow needs them.`,
        icon: "warning",
        showCancelButton: true,
        confirmButtonText: "Remove all",
        cancelButtonText: "Cancel",
        confirmButtonColor: "#d33",
        theme: "dark",
      });
      if (!confirmation.isConfirmed) return;
      this.removeCache(this.cacheEntries);
    },
    removeCache(entries) {
      if (this.cacheLocked) return;
      try {
        removeCacheEntries(entries.map((entry) => entry.id));
      } catch (error) {
        Swal.fire({ title: "Removing cache entries failed", text: error.message, icon: "error", theme: "dark" });
      }
      this.handleRefreshCache();
    },
    cacheStatusInfo(status) {
      return CACHE_STATUS[status] || { label: status, color: "grey", hint: "" };
    },
    formatBytes(bytes) {
      const units = ["B", "KB", "MB", "GB", "TB"];
      let value = bytes;
      let unit = 0;
      while (value >= 1024 && unit < units.length - 1) {
        value /= 1024;
        unit += 1;
      }
      return `${value.toFixed(unit ? 1 : 0)} ${units[unit]}`;
    },
    formatDateTime(isoDate) {
      return isoDate ? new Date(isoDate).toLocaleString() : "";
    },
    formatRelativeDate(isoDate) {
      if (!isoDate) return "";
      const date = new Date(isoDate);
      const startOfToday = new Date().setHours(0, 0, 0, 0);
      const daysAgo = Math.ceil((startOfToday - date.getTime()) / DAY_MS);
      if (daysAgo <= 0) return "Today";
      if (daysAgo === 1) return "Yesterday";
      if (daysAgo < 30) return `${daysAgo} days ago`;
      return date.toLocaleDateString();
    },
    async handleEngineChange(engine) {
      if (!engine || engine === this.selectedEngine) {
        return;
      }
      await this.activateEngineQuietly(engine);
    },
    async handleStartEngine() {
      await this.activateEngineQuietly(this.selectedEngine);
    },
    async handleDetectEngines() {
      try {
        await this.$store.dispatch("bootstrapContainerRuntime");
      } catch (error) {
        console.error("Container engine detection failed:", error);
      }
    },
    // Failures are shown by ContainerEngineDialog (error + Retry), so no extra popup here.
    // Never changes the saved default; only the "Use as default" checkbox does.
    async activateEngineQuietly(engine) {
      try {
        await this.$store.dispatch("activateEngine", { engine, start: true });
      } catch (error) {
        console.error("Engine start failed:", error);
      }
    },
    handleDefaultChange(value) {
      const checked = typeof value === "boolean" ? value : !this.useAsDefault;
      this.$store.dispatch("setDefaultRuntime", checked ? this.selectedEngine : "");
    },
    handleApplyResources() {
      // Win/mac: restart Docker Desktop or Podman machine so CPU/RAM take effect.
      if (this.resourceEngine === "podman") {
        this.restartPodmanMachine();
        return;
      }
      if (this.$store.state.OStype === "Darwin") {
        this.restartDockerMacOS();
        return;
      }
      if (this.$store.state.OStype === "Windows_NT") {
        this.restartDockerWin();
      }
    },
    updateDockerSettings(memoryMiB, processors) {
      // macOS Docker Desktop only. Windows WSL2 CPU/RAM are in ~/.wslconfig.
      const settingsPath = this.dockerSettingsPath;
      if (!settingsPath) {
        console.error("Docker Desktop settings file was not found.");
        return;
      }
      fs.readFile(settingsPath, "utf-8", (err, data) => {
        if (err) {
          console.error(`error reading file ${settingsPath}`);
          return;
        }
        let settingsJSON = JSON.parse(data);
        const cpuKeys = ["Cpus", "cpus", "CpuCount", "cpuCount"].filter((key) =>
          Object.prototype.hasOwnProperty.call(settingsJSON, key)
        );
        const memKeys = ["MemoryMiB", "memoryMiB"].filter((key) =>
          Object.prototype.hasOwnProperty.call(settingsJSON, key)
        );
        (cpuKeys.length ? cpuKeys : ["Cpus"]).forEach((key) => {
          settingsJSON[key] = processors;
        });
        (memKeys.length ? memKeys : ["MemoryMiB"]).forEach((key) => {
          settingsJSON[key] = memoryMiB;
        });
        const updatedSettings = JSON.stringify(settingsJSON, null, 2);
        fs.writeFile(settingsPath, updatedSettings, "utf-8", (writeErr) => {
          if (writeErr) {
            console.error(`Error writing file ${settingsPath}`);
            return;
          }
          console.log(
            `Settings file updated: CPUs: ${processors}, RAM: ${memoryMiB}`
          );
        });
      });
    },
    async restartDockerWin() {
      try {
        // Step 0: Warning prompt
        const confirmation = await Swal.fire({
          title: 'Warning: Docker and WSL Restart Required',
          html: `
            <div style="text-align: left; margin: 15px 0;">
              <p><strong>This operation will:</strong></p>
              <ul style="margin: 10px 0; padding-left: 20px;">
                <li>Stop all running Docker containers</li>
                <li>Shut down all WSL distributions</li>
                <li>Restart Docker Desktop</li>
                <li>Apply new resource settings: <strong>${this.memtotal}GB RAM, ${this.ncpu} CPUs</strong></li>
              </ul>
              <p><strong>⚠️ Any unsaved work in WSL or Docker containers will be lost.</strong></p>
            </div>
          `,
          icon: 'warning',
          showCancelButton: true,
          confirmButtonText: 'Continue',
          cancelButtonText: 'Cancel',
          confirmButtonColor: '#d33',
          cancelButtonColor: '#3085d6',
          theme: 'dark'
        });

        if (!confirmation.isConfirmed) {
          return; // User cancelled
        }

        // WSL2 VM limits live in ~/.wslconfig, not Docker Desktop's JSON.
        writeWslResourceConfig({
          memoryGb: Math.max(1, Number(this.memtotal) || 1),
          processors: Math.max(1, Number(this.ncpu) || 1),
        });
        
        const progressDialog = Swal.fire({
          title: 'Applying WSL Resource Changes',
          html: 'Step 1/5: Stopping Docker Desktop...',
          allowOutsideClick: false,
          showConfirmButton: false,
          theme: 'dark',
          didOpen: () => Swal.showLoading()
        });

        // Step 1: Stop Docker Desktop processes (best-effort).
        // Stop-Process still exits 1 when a name is missing or access is
        // denied, even with -ErrorAction SilentlyContinue — do not abort.
        await this.execPowerShellCommand(
          "Get-Process | Where-Object { $_.ProcessName -match '^(Docker Desktop|com\\.docker\\.)' } | Stop-Process -Force -ErrorAction SilentlyContinue; exit 0",
          { ignoreError: true }
        );
        
        // Step 2: Shutdown WSL
        progressDialog.update({ html: 'Step 2/5: Shutting down WSL...' });
        await this.execPowerShellCommand('wsl --shutdown');
        
        // Step 3: Wait for complete shutdown (8 second rule with countdown)
        for (let i = 8; i > 0; i--) {
          progressDialog.update({ 
            html: `Step 3/5: Waiting for WSL to stop completely...<br><strong style="font-size: 1.2em; color: #1DE9B6;">${i} seconds remaining</strong>` 
          });
          await this.sleep(1000);
        }
        
        // Step 4: Start Docker Desktop (which will start WSL with new config)
        progressDialog.update({ html: 'Step 4/5: Starting Docker Desktop...' });
        const dockerExe = dockerDesktopWindowsExe();
        if (!dockerExe) {
          throw new Error(
            "Docker Desktop was not found. Install Docker Desktop, then retry."
          );
        }
        await this.execPowerShellCommand(`Start-Process '${dockerExe.replace(/'/g, "''")}'`);
        
        // Step 5: Wait for Docker to initialize
        progressDialog.update({ html: 'Step 5/5: Waiting for Docker to initialize...' });
        await this.sleep(5000);
        
        await progressDialog.close();
        await Swal.fire({
          title: 'Resources Updated Successfully',
          text: `Docker Desktop restarted with ${this.memtotal}GB RAM, ${this.ncpu} CPUs`,
          icon: 'success',
          theme: 'dark',
          timer: 4000
        });
        
      } catch (error) {
        await Swal.fire({
          title: 'Update Failed',
          text: error.message,
          icon: 'error',
          theme: 'dark'
        });
      }
    },
    restartDockerMacOS() {
      console.log("starting edit", Math.round(this.memory * 1024));
      this.updateDockerSettings(Math.round(this.memtotal * 1024), this.ncpu);
      const commands = [
        "osascript -e 'quit app \"Docker\"'",
        "pkill -f Docker",
        "sleep 3 && open --background -a Docker",
      ];

      commands.forEach((command) => {
        exec(command, (error, stdout, stderr) => {
          if (error) {
            console.error(`Error executing command: ${command}`, error);
            return;
          }
          console.log(`Command executed: ${command}`);
          console.log(`stdout: ${stdout}`);
          console.log(`stderr: ${stderr}`);
        });
      });
    },
    restartDockerLinux() {
      console.log("in development");
    },
    async restartPodmanMachine() {
      const isWindows = this.$store.state.OStype === "Windows_NT";
      try {
        const confirmation = await Swal.fire({
          title: "Warning: Podman Machine Restart Required",
          html: `
            <div style="text-align: left; margin: 15px 0;">
              <p><strong>This operation will:</strong></p>
              <ul style="margin: 10px 0; padding-left: 20px;">
                <li>Stop the Podman machine</li>
                ${
                  isWindows
                    ? `<li>Update <code>~/.wslconfig</code> with <strong>${this.memtotal}GB RAM, ${this.ncpu} CPUs</strong> (WSL-wide; Podman on Windows cannot set per-machine CPU/RAM)</li>
                <li>Shut down WSL so the new limits apply</li>`
                    : `<li>Apply new resource settings: <strong>${this.memtotal}GB RAM, ${this.ncpu} CPUs</strong></li>`
                }
                <li>Start the Podman machine again</li>
              </ul>
              <p><strong>Any running containers will be stopped.</strong></p>
            </div>
          `,
          icon: "warning",
          showCancelButton: true,
          confirmButtonText: "Continue",
          cancelButtonText: "Cancel",
          confirmButtonColor: "#d33",
          cancelButtonColor: "#3085d6",
          theme: "dark",
        });

        if (!confirmation.isConfirmed) {
          return;
        }

        const progressDialog = Swal.fire({
          title: "Applying Podman Resource Changes",
          html: isWindows
            ? "Updating .wslconfig and restarting WSL/Podman..."
            : "Updating machine CPU and memory...",
          allowOutsideClick: false,
          showConfirmButton: false,
          theme: "dark",
          didOpen: () => Swal.showLoading(),
        });

        const result = await applyPodmanMachineResources({
          cpus: this.ncpu,
          memoryMiB: Math.round(this.memtotal * 1024),
        });

        progressDialog.update({ html: "Waiting for Podman to come back online..." });
        await this.sleep(3000);
        await this.$store.dispatch("probeContainerRuntimes", { force: true });
        await this.$store.dispatch("fetchDockerInfo");
        await progressDialog.close();

        const viaNote =
          result.appliedVia === "wslconfig"
            ? " (via ~/.wslconfig)"
            : "";
        await Swal.fire({
          title: "Resources Updated Successfully",
          text: `Podman machine ${result.machineName} restarted with ${this.memtotal}GB RAM, ${this.ncpu} CPUs${viaNote}`,
          icon: "success",
          theme: "dark",
          timer: 4000,
        });
      } catch (error) {
        await Swal.fire({
          title: "Update Failed",
          text: error.message || String(error),
          icon: "error",
          theme: "dark",
        });
      }
    },

    // Helper method for sequential PowerShell execution.
    // execFile avoids cmd.exe quote-stripping of nested -Command strings.
    execPowerShellCommand(command, { ignoreError = false } = {}) {
      return new Promise((resolve, reject) => {
        execFile(
          "powershell.exe",
          ["-NoProfile", "-NonInteractive", "-Command", command],
          (error, stdout, stderr) => {
            if (error && !ignoreError) {
              reject(new Error(`Command failed: ${command}\n${error.message}`));
              return;
            }
            console.log(`Executed: ${command}`);
            if (stdout) console.log(`Output: ${stdout}`);
            if (error) console.warn(`Ignored PowerShell error: ${error.message}`);
            resolve({ stdout, stderr });
          }
        );
      });
    },

    // Helper for delays
    sleep(ms) {
      return new Promise(resolve => setTimeout(resolve, ms));
    },
  },
  components: {},
};
</script>

<style scoped>
.rm-page {
  display: flex;
  justify-content: center;
  align-items: center;
  min-height: 100vh;
  padding: 24px;
}
/* Frosted glass: blurs and darkens the background photo so text stays readable. */
.rm-card.v-card {
  width: 100%;
  max-width: 900px;
  background-color: rgba(20, 20, 24, 0.72) !important;
  backdrop-filter: blur(14px);
  border: 1px solid rgba(255, 255, 255, 0.08);
  border-radius: 10px !important;
  box-shadow: 0 8px 32px rgba(0, 0, 0, 0.35) !important;
}
.rm-header {
  display: flex;
  align-items: center;
  padding: 14px 18px;
}
.rm-title {
  font-size: 15px;
  font-weight: 600;
  letter-spacing: 0.08em;
}
.rm-status {
  font-size: 13px;
}
.rm-section,
.rm-columns.no-gutters > .rm-section {
  padding: 12px 18px;
}
.rm-col-divider {
  border-right: 1px solid rgba(255, 255, 255, 0.08);
}
.rm-section-title {
  font-size: 11px;
  font-weight: 700;
  letter-spacing: 0.1em;
  text-transform: uppercase;
  opacity: 0.65;
  margin-bottom: 8px;
}
.rm-muted {
  font-size: 12px;
  opacity: 0.7;
}
.rm-break {
  word-break: break-all;
}
.rm-slider-row {
  display: flex;
  align-items: center;
  height: 32px;
}
.rm-slider-label {
  width: 40px;
  font-size: 13px;
}
.rm-slider-value {
  min-width: 76px;
  text-align: right;
  font-size: 13px;
  font-variant-numeric: tabular-nums;
}
.rm-table.v-data-table {
  background-color: transparent !important;
}
.rm-table ::v-deep th {
  font-size: 11.5px !important;
}
.rm-table ::v-deep tbody tr:hover {
  background-color: rgba(255, 255, 255, 0.05) !important;
}
</style>
