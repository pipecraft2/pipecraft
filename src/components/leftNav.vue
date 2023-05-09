<template>
  <v-list dense rounded>
    <v-list-item>
      <v-list-item-content>
        <v-tooltip right>
          <template v-slot:activator="{ on }">
            <div v-on="on">
              <v-btn
                v-if="$store.state.runInfo.active != true"
                style="background-color: #212121"
                block
                :style="
                  $store.state.inputDir != ''
                    ? {
                        borderBottom: 'thin #1DE9B6 solid',
                        borderTop: 'thin white solid',
                        borderLeft: 'thin white solid',
                        borderRight: 'thin white solid',
                      }
                    : {
                        borderBottom: 'thin #E57373 solid',
                        borderTop: 'thin white solid',
                        borderRight: 'thin white solid',
                        borderLeft: 'thin white solid',
                      }
                "
                @click="folderSelect"
              >
                Select workDir
              </v-btn>
              <v-btn
                v-if="$store.state.runInfo.active == true"
                style="background-color: #212121"
                block
              >
                Processing...
              </v-btn>
              <v-progress-linear
                v-if="$store.state.runInfo.active == true"
                indeterminate
                color="#1DE9B6"
              ></v-progress-linear>
            </div>
          </template>
          <div v-if="this.$store.state.inputDir == ''">No files selected!</div>
          <div v-else>
            <div>{{ $store.state.inputDir }}</div>
            <div>
              {{ $store.state.data.dataFormat }} |
              {{ $store.state.data.readType }} |
              {{ $store.state.data.fileFormat }}
            </div>
          </div>
        </v-tooltip>
      </v-list-item-content>
    </v-list-item>
        <v-list-item ripple >
          <v-menu tile dark left >
            <template v-slot:activator="{ on, attrs }">
              <v-list-item-content v-on="on" v-bind="attrs">
                <v-btn

                :style="
                    $route.path.includes('premade')
                                          ? {
                        borderBottom: 'thin #1DE9B6 solid',
                        borderTop: 'thin white solid',
                        borderLeft: 'thin white solid',
                        borderRight: 'thin white solid',
                      }
                    : {
                        borderBottom: 'thin #E57373 solid',
                        borderTop: 'thin white solid',
                        borderRight: 'thin white solid',
                        borderLeft: 'thin white solid',
                      }
                  "
                block
              >
                Select pipeline
              </v-btn>
                
              </v-list-item-content>
            </template>
            <v-list style="padding: 0">
              <v-subheader style="height: 40px">SELECT WORKFLOW</v-subheader>
              <v-divider></v-divider>
              <v-list-item
                ripple
                link
                v-for="(item, index) in $store.state.customWorkflowInfo"
                :key="index"
                @click="push2premade(index)"
              >
                <v-list-item-title>{{
                  index.replace("_", " ")
                }}</v-list-item-title>
              </v-list-item>
            </v-list>
          </v-menu>
        </v-list-item>
    <v-list-item>
      <v-list-item-content v-if="$store.state.runInfo.active == true">
        <v-btn
          block
          style="
            background-color: #212121;
            border: thin #e57373 solid;
            color: red;
          "
          @click="stopWorkflow"
        >
          Stop
        </v-btn>
      </v-list-item-content>
      <v-list-item-content v-else>
        <RunButton />
      </v-list-item-content>
    </v-list-item>
    <v-divider></v-divider>
    <SelectedRoutes />
  </v-list>
</template>

<script>
// var path = require("path");
const Swal = require("sweetalert2");
const slash = require("slash");
const { dialog } = require("@electron/remote");
// const { dialog } = require("electron").remote;

import RunButton from "./RunButton";
import SelectedRoutes from "./SelectedRoutes";
import * as Dockerode from "dockerode";
import os from "os";
var socketPath =
  os.platform() === "win32" ? "//./pipe/docker_engine" : "/var/run/docker.sock";
var dockerode = new Dockerode({ socketPath: socketPath });

export default {
  name: "leftNav",
  components: { RunButton, SelectedRoutes },
  data() {
    return {
      loader: null,
      loading: false,
    };
  },
  watch: {
    loader() {
      const l = this.loader;
      this[l] = !this[l];

      setTimeout(() => (this[l] = false), 9000);

      this.loader = null;
    },
  },
  methods: {
        push2premade(name) {
      if (this.$route.path != `/premade/${name}`) {
        this.$router.push(`/premade/${name}`);
      }
    },
    async stopWorkflow() {
      var container = dockerode.getContainer(
        this.$store.state.runInfo.containerID
      );
      container.remove({ v: true, force: true });
      this.$store.commit("resetRunInfo");
    },
    folderSelect() {
      Swal.mixin({
        input: "select",
        confirmButtonText: "Next &rarr;",
        showCancelButton: true,
        progressSteps: ["1", "2", "3"],
      })
        .queue([
          {
            title: "Sequencing data format",
            inputOptions: {
              demultiplexed: "demultiplexed",
              multiplexed: "multiplexed",
            },
          },
          {
            title: "Sequence files extension",
            inputOptions: {
              Uncompressed: {
                fastq: "*.fastq",
                fasta: "*.fasta",
                fq: "*.fq",
                fa: "*.fa",
                txt: "*.txt",
              },
              Compressed: {
                fastq_gz: "*.fastq.gz",
                fasta_gz: "*.fasta.gz",
                fq_gz: "*.fq.gz",
                fa_gz: "*.fa.gz",
                txt_gz: "*.txt.gz",
              },
            },
          },
          {
            title: "Sequencing read types",
            inputOptions: {
              paired_end: "paired-end",
              single_end: "single-end",
            },
          },
        ])
        .then(async (result) => {
          if (result.value) {
            this.$store.commit("addInputInfo", {
              readType: result.value[2],
              dataFormat: result.value[0],
              fileFormat: result.value[1].replace("_", "."),
            });
            this.$store.commit("toggle_PE_SE_scripts", result.value[2]);
            this.$store.commit("toggle_demux_mux", result.value[0]);
            dialog
              .showOpenDialog({
                title: "Select the folder containing your sequnece files",
                properties: ["openDirectory", "showHiddenFiles"],
              })
              .then((result) => {
                if (typeof result.filePaths[0] !== "undefined") {
                  var correctedPath = slash(result.filePaths[0]);
                  this.$store.commit("addInputDir", correctedPath);
                }
              })
              .catch((err) => {
                console.log(err);
              });
          }
        });
    },
  },
};
</script>

<style scoped>
.v-btn {
  justify-content: center;
}
.custom-loader {
  animation: loader 1s infinite;
  display: flex;
}
@-moz-keyframes loader {
  from {
    transform: rotate(0);
  }
  to {
    transform: rotate(360deg);
  }
}
@-webkit-keyframes loader {
  from {
    transform: rotate(0);
  }
  to {
    transform: rotate(360deg);
  }
}
@-o-keyframes loader {
  from {
    transform: rotate(0);
  }
  to {
    transform: rotate(360deg);
  }
}
@keyframes loader {
  from {
    transform: rotate(0);
  }
  to {
    transform: rotate(360deg);
  }
}
</style>
