<template>
  <div class="container" >
    <v-card width="85%" style="margin: auto; background-color: rgb(0 0 0 / 40%)">
      <v-card-title style="color: black" class="py-10">RESOURCE MANAGER</v-card-title>
      <v-divider></v-divider>
      <v-card-subtitle
        v-if="!isDockerActive"
        style="color: black; display: flex; align-items: center"
      >
        CPU:
        <v-progress-linear
          rounded
          reverse
          color="black"
          buffer-value="0"
          stream
          style="margin-left: 10px; width: 40px"
        ></v-progress-linear>
      </v-card-subtitle>
      <v-card-subtitle v-else style="color: black"
        >CPU: {{ ncpu }}</v-card-subtitle
      >
      <v-slider
        class="px-5"
        v-model="ncpu"
        color="black"
        track-color="black"
        step="1"
        ticks="always"
        tick-size="5"
        min="1"
        :max="processors"
        :tick-labels="tickLabelsCPU"
        height="100"
        style="color: black"
      ></v-slider>
      <v-divider></v-divider>
      <v-card-subtitle
        v-if="!isDockerActive"
        style="color: black; display: flex; align-items: center"
      >
        RAM:
        <v-progress-linear
          rounded
          reverse
          color="black"
          buffer-value="0"
          stream
          style="margin-left: 10px; width: 40px"
        ></v-progress-linear>
      </v-card-subtitle>
      <v-card-subtitle v-else style="color: black">
        RAM:
        {{ memtotal }}
        GB</v-card-subtitle
      >

      <v-slider
        class="px-5"
        color="black"
        v-model="memtotal"
        track-color="black"
        step="1"
        ticks="always"
        tick-size="5"
        min="1"
        :max="memory"
        :tick-labels="tickLabelsMEM"
        height="100"
        style="color: black"
      ></v-slider>
      <div style="display: flex; justify-content: center">
        <v-btn v-if="this.$store.state.OStype != 'Linux'"
          @click="
            $store.state.OStype === 'Darwin'
              ? restartDockerMacOS()
              : $store.state.OStype === 'Linux'
              ? restartDockerLinux()
              : $store.state.OStype === 'Windows_NT'
              ? restartDockerWin()
              : null
          "
          style="margin: auto"
          class="ma-5"
          outlined
          color="black"
        >
          Apply & restart docker
        </v-btn>
      </div>
    </v-card>
  </div>
</template>

<script>
const { exec } = require("child_process");
const fs = require("fs");
const path = require("path");
import os from "os";
import { mapState, mapGetters } from "vuex";
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
      dockerSettingsPath: state => state.systemSpecs.dockerSettings
    }),
    ...mapGetters(['isDockerActive']),
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

        // Step 1: Stop Docker Desktop processes
        await this.execPowerShellCommand("Stop-Process -Name 'Docker Desktop' -Force -ErrorAction SilentlyContinue");
        await this.execPowerShellCommand("Stop-Process -Name 'com.docker.backend' -Force -ErrorAction SilentlyContinue");
        
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

    // Helper method for sequential PowerShell execution
    execPowerShellCommand(command) {
      return new Promise((resolve, reject) => {
        exec(`powershell.exe -Command "${command}"`, (error, stdout, stderr) => {
          if (error) {
            reject(new Error(`Command failed: ${command}\n${error.message}`));
          } else {
            console.log(`✅ Executed: ${command}`);
            if (stdout) console.log(`Output: ${stdout}`);
            resolve({ stdout, stderr });
          }
        });
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
