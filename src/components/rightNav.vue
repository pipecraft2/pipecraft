<template>
  <v-list dense>
    <v-list-item>
      <v-tooltip left nudge-left="10">
        <template v-slot:activator="{ on }">
          <v-list-item-content v-on="on">
            <v-icon :size="50" :color="dockerActive"> mdi-docker </v-icon>
          </v-list-item-content>
        </template>
        <span v-if="dockerActive === '#1DE9B6'">docker desktop is running</span>
        <span v-else-if="dockerActive === '#FF7043'"
          >docker desktop stopped</span
        >
      </v-tooltip>
    </v-list-item>
    <v-divider class="mt-2 mb-2"></v-divider>
    <v-tooltip left nudge-left="10">
      <template v-slot:activator="{ on }">
        <v-list-item
          :disabled="$store.state.runInfo.active"
          v-on="on"
          class="mt-5"
          ripple
          link
        >
          <v-menu tile dark left nudge-left="15" offset-x>
            <template v-slot:activator="{ on, attrs }">
              <v-list-item-content v-on="on" v-bind="attrs">
                <v-icon
                  large
                  :style="
                    $store.getters.singleStepNames.some((v) =>
                      $route.path.includes(v)
                    )
                      ? { color: '#1DE9B6' }
                      : { color: 'white' }
                  "
                  >mdi-tools</v-icon
                >
              </v-list-item-content>
            </template>
            <v-list style="padding: 0">
              <v-subheader style="height: 40px">QUICK TOOLS</v-subheader>
              <v-divider></v-divider>
              <v-list-item
                ripple
                link
                v-for="(item, index) in $store.state.steps"
                :key="index"
                @click="
                  push2route(item.stepName, index),
                    addStep(item, nrOfSelectedSteps)
                "
              >
                <v-list-item-title>{{ item.stepName }}</v-list-item-title>
              </v-list-item>
            </v-list>
          </v-menu>
        </v-list-item>
      </template>
      <span>Quick Tools. Select and use individual tools</span>
    </v-tooltip>
    <v-divider class="mt-2 mb-2"></v-divider>
    <v-list-item
      class="mt-5"
      ripple
      link
      v-for="item in this.items"
      :key="item.title"
    >
      <v-tooltip left nudge-left="10">
        <template v-slot:activator="{ on }">
          <v-list-item-content v-on="on" @click="item.action">
            <v-icon
              :style="
                `/${item.title}` == $route.path ||
                (item.title == 'Debug' && $store.state.data.debugger == true)
                  ? { color: '#1DE9B6' }
                  : { color: 'white' }
              "
              >{{ item.icon }}</v-icon
            >
          </v-list-item-content>
        </template>
        <span>{{ item.tooltip }}</span>
      </v-tooltip>
    </v-list-item>
  </v-list>
</template>

<script>
import os from "os";
const { dialog } = require("@electron/remote");
const slash = require("slash");
const fs = require("fs");
var Docker = require("dockerode");
var socketPath =
  os.platform() === "win32" ? "//./pipe/docker_engine" : "/var/run/docker.sock";
var docker = new Docker({ socketPath: socketPath });
var JSONfn = require("json-fn");

export default {
  name: "rightNav",
  data() {
    return {
      nrOfSelectedSteps: (state) => state.selectedSteps.length + 1,
      dockerActive: "pending",
      items: [
        {
          title: "Save",
          icon: "mdi-content-save",
          tooltip: "save workflow",
          action: this.saveWorkFlow,
        },
        {
          title: "Load",
          icon: "mdi-folder-open",
          tooltip: "load workflow",
          action: this.loadWorkFlow,
        },
        {
          title: "fastqcANDmultiqc",
          icon: "mdi-beaker-check",
          tooltip: "QualityCheck (FastQC and MultiQC)",
          action: this.push2qc,
        },
        {
          title: "ExpertMode",
          icon: "mdi-puzzle-edit",
          tooltip: "Expert mode",
          action: this.push2expert,
        },
        {
          title: "Debug",
          icon: "mdi-bug",
          tooltip:
            "Green = debugging mode ON; PipeCraft will retain all temporary files",
          action: this.debug,
        },
      ],
    };
  },
  created() {
    var self = this;
    setInterval(async function () {
      self.dockerActive = await docker
        .version()
        .then(() => {
          if (self.dockerActive != "#1DE9B6") {
            self.$store.commit("updateDockerStatus", "running");
          }
          return "#1DE9B6";
        })
        .catch(() => {
          self.$store.commit("updateDockerStatus", "stopped");
          return "#FF7043";
        });
    }, 1000);
  },
  methods: {
    addStep(item) {
      this.$store.commit("addStep", {
        step: item,
      });
    },
    saveWorkFlow() {
      dialog
        .showSaveDialog({
          title: "Save current configuration",
          filters: [{ name: "JSON", extensions: ["JSON"] }],
        })
        .then((result) => {
          if (result.canceled !== true) {
            var conf = [];
            let confJson;
            let configSavePath = slash(result.filePath);
            if (this.$route.params.workflowName) {
              conf.push(this.$store.state[this.$route.params.workflowName]);
              conf.push(this.$route.params.workflowName);
              confJson = JSONfn.stringify(conf);
            } else {
              confJson = JSONfn.stringify(this.$store.state.selectedSteps);
            }
            fs.writeFileSync(configSavePath, confJson);
          }
        })
        .catch((err) => {
          console.log(err);
        });
    },
    loadWorkFlow() {
      dialog
        .showOpenDialog({
          title: "Select a previous configuration",
          filters: [{ name: "JSON", extensions: ["JSON"] }],
        })
        .then((result) => {
          if (result.canceled !== true) {
            let configLoadPath = slash(result.filePaths[0]);
            let configJSON = fs.readFileSync(configLoadPath);
            let configObj = JSONfn.parse(configJSON);
            if (Object.keys(this.$store.state).includes(configObj[1])) {
              this.$store.commit("loadCustomWorkflow", configObj);
              this.$router.push(`/premade/${configObj[1]}`);
            } else {
              this.$store.commit("loadWorkflow", configObj);
            }
          }
        });
    },
    debug() {
      this.$store.commit("toggleDebugger");
    },
    push2premade(name) {
      if (this.$route.path != `/premade/${name}`) {
        this.$router.push(`/premade/${name}`);
      }
    },
    push2qc() {
      if (this.$route.path != "/fastqcANDmultiqc") {
        this.$router.push("/fastqcANDmultiqc");
      }
    },
    push2expert() {
      if (this.$route.path != "/ExpertMode") {
        this.$router.push("/ExpertMode");
      }
    },
    push2route(stepName) {
      let route = `/step/${stepName}/0`;
      if (this.$route.path != route) {
        this.$router.push(route);
      }
    },
  },
};
</script>

<style scoped>
.column {
  float: left;
  width: 33.33%;
  padding: 5px;
}

/* Clear floats after image containers */
.row::after {
  content: "";
  clear: both;
  display: table;
}
.v-icon.v-icon::after {
  transform: scale(1);
}
.material icons.primary header material icon first menu {
  margin-left: -2 px;
}
</style>
