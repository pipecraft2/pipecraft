<template>
  <div class="container" >
    <v-card width="85%" style="margin: auto; background-color: rgb(0 0 0 / 40%)">
      <v-card-title style="color: white" class="py-10">RESOURCE MANAGER</v-card-title>
      <v-divider></v-divider>
      <v-card-subtitle style="color: white">Container runtime</v-card-subtitle>
      <div class="px-5 pb-4">
        <!-- Docker | Podman switch when both are installed. Auto is gone. -->
        <v-btn-toggle
          v-if="bothEnginesDetected"
          :value="selectedEngine"
          :mandatory="Boolean(selectedEngine)"
          dense
          color="teal accent-3"
          background-color="transparent"
          :disabled="runtimeStarting"
          @change="handleEngineChange"
        >
          <v-btn value="docker" small dark :disabled="!availableRuntimes.docker">Docker</v-btn>
          <v-btn value="podman" small dark :disabled="!availableRuntimes.podman">Podman</v-btn>
        </v-btn-toggle>
        <v-checkbox
          v-if="bothEnginesDetected"
          dark
          hide-details
          class="mt-3"
          :disabled="!selectedEngine || runtimeStarting"
          :input-value="useAsDefault"
          @change="handleDefaultChange"
          label="Use as default"
        ></v-checkbox>
        <div
          v-if="bothEnginesDetected"
          class="mt-1"
          style="color: white; font-size: 13px; opacity: 0.8"
        >
          When checked, this engine starts automatically next time. Uncheck to choose again at launch.
        </div>
        <div class="mt-3" style="color: white; font-size: 14px">
          <div>{{ runtimeStatusLine }}</div>
          <div v-if="socketPath" style="opacity: 0.8; word-break: break-all">
            Socket: {{ socketPath }}
          </div>
          <div style="opacity: 0.8">
            Detected: {{ detectedRuntimesText }}
          </div>
          <div v-if="!showApplyRestart" class="mt-2" style="opacity: 0.8">
            CPU and RAM limits are applied to each workflow container. No engine restart is required on Linux.
          </div>
        </div>
      </div>
      <v-divider></v-divider>
      <v-card-subtitle
        v-if="!isDockerActive"
        style="color: white; display: flex; align-items: center"
      >
        CPU:
        <v-progress-linear
          rounded
          reverse
          color="white"
          buffer-value="0"
          stream
          style="margin-left: 10px; width: 40px"
        ></v-progress-linear>
      </v-card-subtitle>
      <v-card-subtitle v-else style="color: white"
        >CPU: {{ ncpu }}</v-card-subtitle
      >
      <v-slider
        class="px-5"
        v-model="ncpu"
        color="white"
        track-color="white"
        step="1"
        ticks="always"
        tick-size="5"
        min="1"
        :max="processors"
        :tick-labels="tickLabelsCPU"
        height="100"
        style="color: white"
      ></v-slider>
      <v-divider></v-divider>
      <v-card-subtitle
        v-if="!isDockerActive"
        style="color: white; display: flex; align-items: center"
      >
        RAM:
        <v-progress-linear
          rounded
          reverse
          color="white"
          buffer-value="0"
          stream
          style="margin-left: 10px; width: 40px"
        ></v-progress-linear>
      </v-card-subtitle>
      <v-card-subtitle v-else style="color: white">
        RAM:
        {{ memtotal }}
        GB</v-card-subtitle
      >

      <v-slider
        class="px-5"
        color="white"
        v-model="memtotal"
        track-color="white"
        step="1"
        ticks="always"
        tick-size="5"
        min="1"
        :max="memory"
        :tick-labels="tickLabelsMEM"
        height="100"
        style="color: white"
      ></v-slider>
      <div style="display: flex; justify-content: center">
        <v-btn v-if="showApplyRestart"
          @click="handleApplyResources"
          style="margin: auto"
          class="ma-5"
          outlined
          color="white"
        >
          {{ applyButtonText }}
        </v-btn>
      </div>
    </v-card>
  </div>
</template>

<script>
const { exec, execFile } = require("child_process");
const fs = require("fs");
const path = require("path");
import os from "os";
import { mapState, mapGetters } from "vuex";
import { applyPodmanMachineResources } from "../utils/containerRuntime";
const CPU = os.cpus().length;
const MEM = Number((os.totalmem() / 1024 ** 3).toFixed(0));
const homeDir = require("os").homedir();
// Path to the .wslconfig file
const wslConfigPath = path.join(homeDir, ".wslconfig");
const Swal = require("sweetalert2");

const createNumberList = (start, end) =>
  Array.from({ length: end - start + 1 }, (_, i) => start + i);

function modTickLabels(arr) {
  const modValue = arr.length >= 48 ? 8 : arr.length > 12 ? 4 : 2;
  return arr.map((value, index) =>
    index === 0 || value % modValue === 0 ? value : ""
  );
}

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
    applyButtonText() {
      return this.resourceEngine === "podman"
        ? "Apply & restart Podman machine"
        : "Apply & restart Docker";
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
  },
  data() {
    return {
      processors: CPU,
      memory: MEM,
      tickLabelsCPU: modTickLabels(createNumberList(1, CPU)),
      tickLabelsMEM: modTickLabels(createNumberList(1, MEM)),
    };
  },
  methods: {
    async handleEngineChange(engine) {
      if (!engine || engine === this.selectedEngine) {
        return;
      }
      try {
        await this.$store.dispatch("activateEngine", {
          engine,
          persistDefault: this.useAsDefault,
          start: true,
        });
      } catch (error) {
        await Swal.fire({
          title: "Engine start failed",
          text: error.message || String(error),
          icon: "error",
          theme: "dark",
        });
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
        updateWslConfig(memory, processors) {
      // Check if .wslconfig file exists
      if (fs.existsSync(wslConfigPath)) {
        // File exists, read and update it
        fs.readFile(wslConfigPath, "utf8", (err, data) => {
          if (err) {
            console.error("Error reading .wslconfig file:", err);
            return;
          }

          // Update memory and processors values
          let updatedConfig = data.replace(/memory=\d+GB/i, `memory=${memory}`);
          updatedConfig = updatedConfig.replace(/processors=\d+/i, `processors=${processors}`);

          // Write the updated content back to the file
          fs.writeFile(wslConfigPath, updatedConfig, "utf8", (err) => {
            if (err) {
              console.error("Error writing to .wslconfig file:", err);
              return;
            }
            console.log(`.wslconfig updated: memory=${memory}, processors=${processors}`);
          });
        });
      } else {
        // File doesn't exist, create it
        const defaultConfig = `[wsl2]
memory=${memory}
processors=${processors}
`;
        fs.writeFile(wslConfigPath, defaultConfig, "utf8", (err) => {
          if (err) {
            console.error("Error creating .wslconfig file:", err);
            return;
          }
          console.log(`.wslconfig created: memory=${memory}, processors=${processors}`);
        });
      }
    },
    updateDockerSettings(memory, processors) {
      fs.readFile(this.dockerSettingsPath, "utf-8", (err, data) => {
        if (err) {
          console.error(`error reading file ${this.dockerSettingsPath}`);
          return;
        }
        let settingsJSON = JSON.parse(data);
        settingsJSON.Cpus = processors;
        settingsJSON.MemoryMiB = memory;
        const updatedSettings = JSON.stringify(settingsJSON, null, 2);
        fs.writeFile(this.dockerSettingsPath, updatedSettings, "utf-8", (err) => {
          if (err) {
            console.error(`Error writing file ${this.dockerSettingsPath}`);
            return;
          }
          console.log(
            `Settings file updated: CPUs: ${processors}, RAM: ${memory}`
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

        // Update WSL config first
        this.updateWslConfig(`${this.memtotal}GB`, this.ncpu);
        
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
        await this.execPowerShellCommand("Start-Process 'C:\\Program Files\\Docker\\Docker\\Docker Desktop.exe'");
        
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
.container {
  display: flex;
  justify-content: center; /* Center horizontally */
  align-items: center; /* Center vertically */
  height: 100vh; /* Full height of the viewport for demonstration */
}
</style>
