<template>
  <v-container pa-10 fluid>
    <div
      style="
        overflow-y: scroll !important;
        min-width: 800px !important;
        min-height: 40% !important;
        height: 60% !important;
        background-color: #454442;
      "
      id="terminal"
    ></div>
    <div
      style="
        min-width: 800px !important;
        background-color: #212121;
        color: white;
      "
    >
      <!-- <div style="padding: 6px">All bioinformatic tools used by pipecraft are available as docker images and can be used via dockers CLI.</div> -->
      <!-- <div style="padding: 6px">Check them out on dockerhub: <a href="https://hub.docker.com/u/pipecraft" target="_blank">https://hub.docker.com/u/pipecraft</a></div>
    <div style="padding: 6px">Some examples for using pipecrafts images with docker CLI are shown below</div>
    <v-divider style="background-color:white"></v-divider>
    <div style="padding: 6px">docker pull pipecraft/dada2:1.20</div>
    <div style="padding: 6px">docker run --interactive --tty -v C:\Users\Name\MyFiles\:/MyFilesOnContainer pipecraft/dada2:1.20</div>
    <v-divider style="background-color:white"></v-divider>
    <div style="padding: 6px">Check out the <a href="https://docs.docker.com/engine/reference/commandline/run/" target="_blank">Docker documentation</a> for more info</div> -->
      <v-card>
        <v-toolbar flat color="grey" dark>
          <v-toolbar-title
            >Using Pipecrafts assets with Docker CLI</v-toolbar-title
          >
        </v-toolbar>
        <v-tabs color="rgb(29, 233, 182)" vertical>
          <v-tab> Managing assets </v-tab>
          <v-tab> Starting a container </v-tab>
          <v-tab> Examples </v-tab>

          <v-tab-item style="background: grey">
            <v-container>
              <v-row justify="space-around">
                <v-card
                  style="
                    min-height: 300px;
                    border-left: thin grey solid;
                    border-radius: 0;
                  "
                >
                  <v-card-subtitle
                    ><div class="font-weight-bold ml-8 mb-2">
                      Bioinformatic tools used by Pipecraft are stored on
                      Dockerhub:
                      <a
                        href="https://hub.docker.com/u/pipecraft"
                        target="_blank"
                        >hub.docker.com/u/pipecraft</a
                      >
                      as Docker images. These images can be used to launch any
                      tool manually with the Docker CLI and even further
                      customize your operations.
                    </div></v-card-subtitle
                  >
                  <v-divider></v-divider>
                  <v-card-text class="pa-0 ma-0">
                    <v-row no-gutters>
                      <v-col
                        v-for="(item, i) in assets"
                        offset="1"
                        cols="8"
                        :key="i"
                      >
                        <div class="font-weight-light">
                          {{ item.info }}
                        </div>
                        <div class="pb-1 font-italic font-weight-medium">
                          {{ item.cmd }}
                          <v-icon
                            @click="copy2clip(item.cmd)"
                            small
                            class="ml-2"
                            >mdi-content-copy</v-icon
                          >
                        </div>
                      </v-col>
                    </v-row>
                  </v-card-text>
                </v-card>
              </v-row>
            </v-container>
          </v-tab-item>
          <v-tab-item>
            <v-container>
              <v-row justify="space-around">
                <v-card
                  style="
                    min-height: 300px;
                    border-left: thin grey solid;
                    border-radius: 0;
                  "
                >
                  <v-card-subtitle
                    ><div class="font-weight-bold ml-8 mb-2">
                      Using Pipecrafts Docker images requires them to be ran as
                      containers which can be achieved with the run command
                      <a
                        href="https://docs.docker.com/engine/reference/run/"
                        target="_blank"
                        >docs.docker.com/engine/reference/run</a
                      >
                      .
                    </div></v-card-subtitle
                  >
                  <v-divider></v-divider>
                  <v-card-text class="pa-0 ma-0">
                    <v-row no-gutters>
                      <v-col
                        v-for="(item, i) in flags"
                        offset="1"
                        cols="9"
                        :key="i"
                      >
                        <div class="font-weight-light">
                          {{ item.info }}
                        </div>
                        <div class="pb-1 font-italic font-weight-medium">
                          {{ item.cmd }}
                          <v-icon
                            @click="copy2clip(item.cmd)"
                            small
                            class="ml-2"
                            >mdi-content-copy</v-icon
                          >
                        </div>
                      </v-col>
                    </v-row>
                  </v-card-text>
                </v-card>
              </v-row>
            </v-container>
          </v-tab-item>
          <v-tab-item>
            <v-container>
              <v-row justify="space-around">
                <v-card
                  style="
                    min-height: 300px;
                    border-left: thin grey solid;
                    border-radius: 0;
                  "
                >
                  <v-card-subtitle
                    ><div class="font-weight-bold ml-8 mb-2">
                      They are many different ways for interacting with
                      containers the easiest way to get started is to attach
                      your files with the -v flag and to keep open an
                      interactive terminal with the --tty and -i flags.
                    </div></v-card-subtitle
                  >
                  <v-divider></v-divider>
                  <v-card-text class="pa-0 ma-0">
                    <v-row no-gutters>
                      <v-col
                        v-for="(item, i) in examples"
                        offset="1"
                        cols="11"
                        :key="i"
                      >
                        <div class="font-weight-light">
                          {{ item.info }}
                        </div>
                        <div class="pb-1 font-italic font-weight-medium">
                          {{ item.cmd }}
                          <v-icon
                            @click="copy2clip(item.cmd)"
                            small
                            class="ml-2"
                            >mdi-content-copy</v-icon
                          >
                        </div>
                      </v-col>
                    </v-row>
                  </v-card-text>
                </v-card>
              </v-row>
            </v-container></v-tab-item
          >
        </v-tabs>
      </v-card>
    </div>
  </v-container>
</template>

<script>
const pty = require("@electron/remote").require("node-pty-prebuilt-multiarch");
const { clipboard } = require("electron");
import os from "os";
import "xterm/css/xterm.css";
import { Terminal } from "xterm";
import { FitAddon } from "xterm-addon-fit";
var term;
var shell = os.platform() === "win32" ? "powershell.exe" : "bash";
var ptyProc = pty.spawn(shell, [], {
  name: "xterm-color",
});
term = new Terminal({
  theme: { background: "#454442" },
  rows: 35,
  experimentalCharAtlas: "dynamic",
});

const fitAddon = new FitAddon();
term.loadAddon(fitAddon);
term.onData((data) => {
  ptyProc.write(data);
});
ptyProc.on("data", function (data) {
  term.write(data);
});
term.onResize((size) => {
  ptyProc.resize(
    Math.max(size ? size.cols : term.cols, 1),
    Math.max(size ? size.rows : term.rows, 1)
  );
  fitAddon.fit();
});

export default {
  created() {
    window.addEventListener("resize", this.myEventHandler);
  },
  destroyed() {
    window.removeEventListener("resize", this.myEventHandler);
  },
  data() {
    return {
      term: null,
      assets: [
        { info: "Show a list of all Images", cmd: "docker images" },
        { info: "Download an image", cmd: "docker pull REPO/IMAGE:TAG" },
        { info: "Delete an image", cmd: "docker rmi IMAGE" },
        { info: "Show a list of all containers", cmd: "docker ps -a" },
        { info: "Delete a running container", cmd: "docker rm -f CONTAINER" },
        { info: "Delete stopped containers", cmd: "docker container prune" },
      ],
      flags: [
        {
          info: "General form of a Docker run command",
          cmd: "docker run [OPTIONS] REPOSITORY/IMAGE:TAG [COMMAND]",
        },
        {
          info: "Attache a folder from the host system",
          cmd: "--volume=host-src:container-dest",
        },
        {
          info: "Set environment variables in the container",
          cmd: '-e ""foo=bar',
        },
        { info: "Keep STDIN open", cmd: "--interactive" },
        { info: "Allocate a psedu terminal", cmd: "--tty" },
      ],
      examples: [
        {
          info: "Run a dada2 container with a active terminal and a attached host folder",
          cmd: "docker run -i --tty -v users/Tom/Files/:/Files pipecraft/dada2:1.20",
        },
        {
          info: "Run a vsearch container with an attached host folder and exectue a script from it",
          cmd: "docker run -v users/Tom/Files/:/Files pipecraft/vsearch:2.18 ./Files/script.sh",
        },
      ],
    };
  },
  mounted() {
    term.open(document.getElementById("terminal"));
    fitAddon.fit();
  },
  beforeDestroy() {},
  name: "ExperMode",
  components: {
    // SelectedRoutes,
  },
  methods: {
    myEventHandler() {
      fitAddon.fit();
    },
    copy2clip(cmd) {
      clipboard.writeText(cmd);
      console.log(cmd);
    },
  },
};
</script>

<style scoped></style>
