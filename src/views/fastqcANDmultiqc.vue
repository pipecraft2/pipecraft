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
        <v-img class="white--text align-end" src="../assets/fastqc_logo.png">
        </v-img>
      </div>
      <div class="spacer" ></div>
      <div class="loader-column" v-if="$store.state.pullLoader.active">
        <div class="image-loader-container">
          <v-progress-circular
            :size="100"
            :width="7"
            color="orange"
            indeterminate
          >
            <span class="loading-text">Pulling<br>image</span>
          </v-progress-circular>
        </div>
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
      <v-tooltip bottom :disabled="isDockerActive && folderPath != ''">
        <template v-slot:activator="{ on }">
          <div v-on="on">
            <v-btn
              :disabled="!isDockerActive || reportLoading || folderPath == ''"
              @click="fastQualityCheck()"
              color="orange"
              text
            >
              Create Report
            </v-btn>
          </div>
        </template>
        <div v-if="!isDockerActive">Failed to find Docker</div>
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
import { mapGetters } from 'vuex';
const shell = require("electron").shell;
const streams = require("memory-streams");
var stdout = new streams.WritableStream();
var stderr = new streams.WritableStream();

export default {
  name: "fastqcANDmultiqc",
  computed: {
    ...mapGetters(['isDockerActive'])
  },
  data() {
    return this.$store.state.Qcheck;
  },
  methods: {
    async folderSelect() {
      await this.$store.dispatch('setWorkingDir', 'fastqcANDmultiqc');
    },
    async fastQualityCheck() {
      this.$store.state.Qcheck.reportReady = false;
      this.$store.state.Qcheck.reportLoading = true;
      console.log("starting fastqc");
      await this.$store.dispatch('imageCheck', "staphb/fastqc:0.11.9");
      await this.$store.dispatch('imageCheck', "ewels/multiqc:1.10");
      let result = await this.$docker
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
      let result2 = await this.$docker
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

<style scoped>
.image-loader-container {
  display: flex;
  justify-content: center;
  align-items: center;
  height: 100%;
}

.loading-text {
  font-size: 14px;
  text-align: center;
  color: white;
  line-height: 1.2;
}

.row {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  position: relative;
  width: 100%;
}

.column {
  flex: 0 0 auto;
  padding: 0 10px;
}

.spacer {
  flex-grow: 1;
}

.loader-column {
  padding-right: 30px; /* Add some space from the right edge */
}
</style>
