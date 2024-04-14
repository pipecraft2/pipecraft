<template>
  <v-card
    class="mx-auto"
    max-width="90%"
    style="margin-top: 200px; background-color: grey; color: white"
  >
    <div class="row" style="padding-left: 25px; padding-top: 25px">
      <div class="column">
        <v-img class="white--text align-end" src="../assets/MultiQC_logo.png">
        </v-img>
      </div>
      <div class="column">
        <v-img class="white--text align-end" src="../assets/fastqc_icon.png">
        </v-img>
      </div>
    </div>

    <v-card-title>FastQC and MultiQC</v-card-title>
    <v-card-subtitle style="color: white" class="pb-0">
      Check out their documentation for more info
    </v-card-subtitle>
    <v-divider class="mt-1"></v-divider>

    <v-card-text class="text--primary">
      <div>
        <a style="color: white" :href="'https://multiqc.info/'" target="_blank"
          >multiqc.info</a
        >
      </div>
      <div>
        <a
          style="color: white"
          :href="'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/'"
          target="_blank"
          >bioinformatics.babraham.ac.uk/projects/fastqc</a
        >
      </div>
    </v-card-text>
    <v-divider class="mt-1"></v-divider>
    <v-card-actions>
      <v-btn
        :disabled="reportLoading"
        @click="folderSelect()"
        color="orange"
        text
      >
        Select Folder
      </v-btn>
      <v-tooltip bottom :disabled="dockerActive && folderPath != ''">
        <template v-slot:activator="{ on }">
          <div v-on="on">
            <v-btn
              :disabled="!dockerActive || reportLoading || folderPath == ''"
              @click="fastQualityCheck()"
              color="orange"
              text
            >
              Create Report
            </v-btn>
          </div>
        </template>
        <div v-if="!dockerActive">Failed to find Docker</div>
        <div v-if="folderPath == ''">No folder selected</div>
      </v-tooltip>
      <v-tooltip right :disabled="reportReady">
        <template v-slot:activator="{ on }">
          <div v-on="on">
            <v-btn
              @click="openReport()"
              color="orange"
              text
              :disabled="!reportReady"
              :loading="reportLoading"
            >
              View Report
              <template v-slot:loader>
                <span>Loading...</span>
              </template>
            </v-btn>
          </div>
        </template>
        <div>No reports generated</div>
      </v-tooltip>
    </v-card-actions>
    <v-progress-linear
      :active="reportLoading"
      color="orange"
      indeterminate
      reverse
    ></v-progress-linear>
  </v-card>
</template>

<script>
import { pullImageAsync } from "dockerode-utils";
import { imageExists } from "dockerode-utils";
import * as Dockerode from "dockerode";
import os from "os";
const shell = require("electron").shell;
const { dialog } = require("@electron/remote");
const slash = require("slash");
const streams = require("memory-streams");
const Swal = require("sweetalert2");
var stdout = new streams.WritableStream();
var stderr = new streams.WritableStream();
var socketPath =
  os.platform() === "win32" ? "//./pipe/docker_engine" : "/var/run/docker.sock";
var dockerode = new Dockerode({ socketPath: socketPath });

export default {
  name: "qualityPlots",
  data() {
    return this.$store.state.Qcheck;
  },
  created() {
    var self = this;
    setInterval(async function () {
      self.dockerActive = await dockerode
        .version()
        .then(() => {
          return true;
        })
        .catch(() => {
          return false;
        });
    }, 1000);
  },
  methods: {
    folderSelect() {
      this.$store.state.Qcheck.reportReady = false;
      Swal.mixin({
        input: "select",
        confirmButtonText: "Next &rarr;",
        showCancelButton: true,
      })
        .queue([
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
        ])
        .then(async (result) => {
          if (result.value) {
            console.log(result.value);
            this.$store.state.Qcheck.fileExtension = result.value[0].replace(
              "_",
              "."
            );
            dialog
              .showOpenDialog({
                title: "Select the folder containing your sequnece files",
                properties: ["openDirectory", "showHiddenFiles"],
              })
              .then((result) => {
                if (typeof result.filePaths[0] !== "undefined") {
                  this.$store.state.Qcheck.folderPath = slash(
                    result.filePaths[0]
                  );
                }
              })
              .catch((err) => {
                console.log(err);
              });
          }
        });
    },
    async fastQualityCheck() {
      this.$store.state.Qcheck.reportReady = false;
      this.$store.state.Qcheck.reportLoading = true;
      console.log("starting fastqc");
      let gotImg = await imageExists(dockerode, "staphb/fastqc:0.11.9");
      if (gotImg === false) {
        console.log(`Pulling image staphb/fastqc:0.11.9`);
        let output = await pullImageAsync(dockerode, "staphb/fastqc:0.11.9");
        console.log(output);
        console.log(`Pull complete`);
      }
      gotImg = await imageExists(dockerode, "ewels/multiqc:1.10");
      if (gotImg === false) {
        console.log(`Pulling image ewels/multiqc:1.10`);
        let output = await pullImageAsync(dockerode, "ewels/multiqc:1.10");
        console.log(output);
        console.log(`Pull complete`);
      }
      let result = await dockerode
        .run(
          "staphb/fastqc:0.11.9",
          [
            "sh",
            "-c",
            `mkdir quality_check | fastqc --outdir quality_check *$format`,
          ],
          [stdout, stderr],
          {
            Tty: false,
            WorkingDir: "/input",
            HostConfig: {
              Binds: [`${this.$store.state.Qcheck.folderPath}:/input`],
            },
            Env: [`format=${this.$store.state.Qcheck.fileExtension}`],
          }
        )
        .then(async ([res, container]) => {
          console.log(stdout.toString());
          let resObj = { statusCode: res.StatusCode };
          container.remove();
          if (res.StatusCode === 0) {
            resObj.log = stdout.toString();
            return resObj;
          } else {
            resObj.log = stderr.toString();
            return resObj;
          }
        })
        .catch((err) => {
          let resObj = {};
          resObj.statusCode = err.statusCode;
          resObj.log = err.json.message;
          return resObj;
        });
      console.log(result);
      stdout = new streams.WritableStream();
      stderr = new streams.WritableStream();
      console.log("starting multiqc");
      let result2 = await dockerode
        .run("ewels/multiqc:1.10", [], [stdout, stderr], {
          Tty: false,
          WorkingDir: "/input",
          HostConfig: {
            Binds: [
              `${this.$store.state.Qcheck.folderPath}/quality_check:/input`,
            ],
          },
          Env: [`format=${this.$store.state.Qcheck.fileExtension}`],
        })
        .then(async ([res, container]) => {
          console.log(stdout.toString());
          let resObj = { statusCode: res.StatusCode };
          container.remove();
          if (res.StatusCode === 0) {
            resObj.log = stdout.toString();
            return resObj;
          } else {
            resObj.log = stderr.toString();
            return resObj;
          }
        })
        .catch((err) => {
          let resObj = {};
          resObj.statusCode = err.statusCode;
          resObj.log = err.json.message;
          return resObj;
        });
      console.log(result2);
      this.$store.state.Qcheck.reportReady = true;
      this.$store.state.Qcheck.reportLoading = false;
    },
    openReport() {
      shell.openExternal(
        `file://${this.$store.state.Qcheck.folderPath}/quality_check/multiqc_report.html`
      );
    },
  },
};
</script>

<style></style>
