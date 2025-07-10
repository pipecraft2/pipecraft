<template>
  <div class="container">
    <v-card width="85%" style="margin: auto">
      <v-card-title class="py-10">RESOURCE MANAGER</v-card-title>
      <v-divider></v-divider>
      <v-card-subtitle
        v-if="this.$store.state.dockerStatus == 'stopped'"
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
      ></v-slider>
      <v-divider></v-divider>
      <v-card-subtitle
        v-if="this.$store.state.dockerStatus == 'stopped'"
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
const CPU = os.cpus().length;
const MEM = Number((os.totalmem() / 1024 ** 3).toFixed(0));
console.log(CPU, MEM);
const homeDir = require("os").homedir();
console.log(homeDir);
// Path to the .wslconfig file
const wslConfigPath = path.join(homeDir, ".wslconfig");
const dockerSettingsPath = [
  `${homeDir}/Library/Group Containers/group.com.docker/settings.json`,
  "/Library/Group Containers/group.com.docker/settings.json",
  `${homeDir}/Library/Containers/com.docker.docker/Data/database/com.docker.driver.amd64-linux/settings.json`,
  `${homeDir}/Library/Containers/com.docker.docker/Data/database/com.docker.driver.amd64-linux/config/daemon.json`,
  `${homeDir}/Library/Containers/com.docker.docker/Data/settings.json`,
  `${homeDir}/Library/Containers/com.docker.docker/Data/com.docker.driver.amd64-linux/settings.json`,
  `${homeDir}/Library/Application Support/Docker Desktop/settings.json`,
  `${homeDir}/Library/Group Containers/group.com.docker/settings-store.json`,
  "/Applications/Docker.app/Contents/Resources/settings.json"
].find(fs.existsSync);
console.log("Docker settings path found:", dockerSettingsPath);

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
      // Read the existing .wslconfig file
      fs.readFile(wslConfigPath, "utf8", (err, data) => {
        if (err) {
          console.error("Error reading .wslconfig file:", err);
          return;
        }

        // Use regular expressions to update memory and processors values
        let updatedConfig = data.replace(/memory=\d+GB/i, `memory=${memory}`);
        updatedConfig = updatedConfig.replace(
          /processors=\d+/i,
          `processors=${processors}`
        );

        // Write the updated content back to the file
        fs.writeFile(wslConfigPath, updatedConfig, "utf8", (err) => {
          if (err) {
            console.error("Error writing to .wslconfig file:", err);
            return;
          }
          console.log(
            `.wslconfig updated: memory=${memory}, processors=${processors}`
          );
        });
      });
    },
    updateDockerSettings(memory, processors) {
      fs.readFile(dockerSettingsPath, "utf-8", (err, data) => {
        if (err) {
          console.error(`error reading file ${dockerSettingsPath}`);
          return;
        }
        let settingsJSON = JSON.parse(data);
        settingsJSON.Cpus = processors;
        settingsJSON.MemoryMiB = memory;
        const updatedSettings = JSON.stringify(settingsJSON, null, 2);
        fs.writeFile(dockerSettingsPath, updatedSettings, "utf-8", (err) => {
          if (err) {
            console.error(`Error writing file ${dockerSettingsPath}`);
            return;
          }
          console.log(
            `Settings file updated: CPUs: ${processors}, RAM: ${memory}`
          );
        });
      });
    },
    restartDockerWin() {
      this.updateWslConfig(`${this.memtotal}GB`, this.ncpu);
      const commands = [
        'Stop-Process -Name "Docker*Desktop" -Force',
        'Stop-Process -Name "com.docker.backend" -Force',
        "wsl --shutdown; Start-Sleep -Seconds 3; wsl",
        "Start-Process '\"C:\\Program Files\\Docker\\Docker\\Docker Desktop.exe\"'",
      ];

      commands.forEach((command) => {
        exec(
          `powershell.exe -Command "${command}"`,
          (error, stdout, stderr) => {
            if (error) {
              console.error(`Error executing command: ${command}`, error);
              return;
            }
            console.log(`Command executed: ${command}`);
            console.log(`stdout: ${stdout}`);
            console.log(`stderr: ${stderr}`);
          }
        );
      });
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
