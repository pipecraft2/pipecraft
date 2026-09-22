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
        color="primary"
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
              color="primary"
              text
            >
              Create Report
            </v-btn>
          </div>
        </template>
        <div v-if="!isDockerActive">{{ engineNotFoundMessage }}</div>
        <div v-if="folderPath == ''">No folder selected</div>
      </v-tooltip>
      <v-tooltip right :disabled="reportReady">
        <template v-slot:activator="{ on }">
          <div v-on="on">
            <v-btn
              @click="openReport()"
              color="primary"
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
      color="primary"
      indeterminate
      reverse
    ></v-progress-linear>
  </v-card>
</template>

<script>
import { mapGetters } from "vuex";
import { PassThrough } from "stream";
import Swal from "sweetalert2";
import {
  prepareBindMounts,
  applyEngineHostConfig,
  getContainerUser,
  hostPathsFromBinds,
  reclaimBindMountOwnership,
} from "../utils/containerRuntime";

const { shell } = require("electron");

export default {
  name: "fastqcANDmultiqc",
  computed: {
    ...mapGetters(["isDockerActive", "engineNotFoundMessage"]),
  },
  data() {
    return this.$store.state.Qcheck;
  },
  methods: {
    async folderSelect() {
      await this.$store.dispatch("setWorkingDir", "fastqcANDmultiqc");
    },

    async fastQualityCheck() {
      this.$store.state.Qcheck.reportReady = false;
      this.$store.state.Qcheck.reportLoading = true;

      const folder = this.$store.state.Qcheck.folderPath;
      const format = this.$store.state.Qcheck.fileExtension;
      const userId = this.$store.state.systemSpecs.userId;
      const groupId = this.$store.state.systemSpecs.groupId;

      try {
        const fastqcResult = await this.runContainer({
          imageName: "staphb/fastqc:0.11.9",
          containerName: "PipeCraft_FastQC",
          workingDir: "/input",
          command: [
            "sh",
            "-c",
            "mkdir -p quality_check && fastqc --outdir quality_check *$format",
          ],
          binds: [`${folder}:/input`],
          env: [
            `format=${format}`,
            `HOST_UID=${userId}`,
            `HOST_GID=${groupId}`,
          ],
          userId,
          groupId,
        });

        if (fastqcResult.StatusCode !== 0) {
          await this.showRunError("FastQC failed", fastqcResult);
          return;
        }

        const multiqcResult = await this.runContainer({
          imageName: "ewels/multiqc:1.10",
          containerName: "PipeCraft_MultiQC",
          workingDir: "/input",
          // Omit command so the image ENTRYPOINT/CMD run as usual.
          binds: [`${folder}/quality_check:/input`],
          env: [
            `format=${format}`,
            `HOST_UID=${userId}`,
            `HOST_GID=${groupId}`,
          ],
          userId,
          groupId,
        });

        if (multiqcResult.StatusCode !== 0) {
          await this.showRunError("MultiQC failed", multiqcResult);
          return;
        }

        this.$store.state.Qcheck.reportReady = true;
      } catch (error) {
        console.error("FastQC/MultiQC error:", error);
        await Swal.fire({
          title: "Quality check failed",
          text: this.errorToMessage(error),
          confirmButtonText: "OK",
          theme: "dark",
        });
      } finally {
        this.$store.state.Qcheck.reportLoading = false;
      }
    },

    /**
     * Same create → attach → demux → wait → remove path as Run.vue.
     */
    async runContainer(spec) {
      const {
        imageName,
        containerName,
        command,
        env = [],
        binds,
        workingDir,
        userId,
        groupId,
      } = spec;

      await this.$store.dispatch("imageCheck", imageName);
      await this.$store.dispatch("clearContainerConflicts", containerName);

      const memory = this.$store.state.dockerInfo.MemTotal;
      const nanoCpus = Math.round(
        Number(this.$store.state.dockerInfo.NCPU) * 1e9
      );
      const createConfig = {
        Image: imageName,
        name: containerName,
        Tty: false,
        AttachStdout: true,
        AttachStderr: true,
        Platform: "linux/amd64",
        Env: env,
        HostConfig: applyEngineHostConfig({
          Binds: prepareBindMounts(binds),
          Memory: memory,
          NanoCpus: nanoCpus,
        }),
      };
      const containerUser = getContainerUser(userId, groupId);
      if (containerUser) {
        createConfig.User = containerUser;
      }
      if (workingDir) {
        createConfig.WorkingDir = workingDir;
      }
      if (Array.isArray(command) && command.length > 0) {
        createConfig.Cmd = command;
      }

      let container = null;
      try {
        container = await this.$docker.createContainer(createConfig);

        const attachStream = await container.attach({
          stream: true,
          stdout: true,
          stderr: true,
        });
        const stdoutStream = new PassThrough();
        const stderrStream = new PassThrough();
        container.modem.demuxStream(attachStream, stdoutStream, stderrStream);

        const endStreams = () => {
          try {
            stdoutStream.end();
          } catch (err) {
            console.debug("stdoutStream.end failed:", err && err.message);
          }
          try {
            stderrStream.end();
          } catch (err) {
            console.debug("stderrStream.end failed:", err && err.message);
          }
        };
        attachStream.on("end", endStreams);
        attachStream.on("close", endStreams);

        const logPromise = this.handleDemuxedStreams(stdoutStream, stderrStream);

        await container.start();
        const data = await container.wait();
        endStreams();

        let stdout = "";
        let stderr = "";
        try {
          const res = await this.waitWithTimeout(logPromise, 2000);
          stdout = res.stdout || "";
          stderr = res.stderr || "";
        } catch (err) {
          console.debug("log drain timeout or error:", err && err.message);
        }

        return {
          StatusCode: data.StatusCode,
          stdout,
          stderr,
        };
      } finally {
        if (container) {
          try {
            await container.remove({ v: true, force: true });
          } catch (err) {
            const msg = err && err.message ? err.message : "";
            if (
              !(msg.includes("HTTP code 404") || msg.includes("HTTP code 409"))
            ) {
              console.warn("Non-fatal remove error:", err);
            }
          }
        }
        try {
          await reclaimBindMountOwnership(hostPathsFromBinds(binds));
        } catch (err) {
          console.warn("Ownership reclaim failed:", err && err.message);
        }
      }
    },

    handleDemuxedStreams(stdoutStream, stderrStream) {
      return new Promise((resolve) => {
        let stdout = "";
        let stderr = "";

        stdoutStream.on("data", (data) => {
          const text = data.toString();
          console.log(text);
          stdout += text;
        });
        stderrStream.on("data", (data) => {
          const text = data.toString();
          console.log(text);
          stderr += text;
        });

        let endedStdout = false;
        let endedStderr = false;
        const tryResolve = () => {
          if (endedStdout && endedStderr) {
            resolve({ stdout, stderr });
          }
        };

        stdoutStream.on("end", () => {
          endedStdout = true;
          tryResolve();
        });
        stderrStream.on("end", () => {
          endedStderr = true;
          tryResolve();
        });
      });
    },

    waitWithTimeout(promise, ms) {
      return new Promise((resolve, reject) => {
        const t = setTimeout(() => reject(new Error("timeout")), ms);
        promise
          .then((v) => {
            clearTimeout(t);
            resolve(v);
          })
          .catch((e) => {
            clearTimeout(t);
            reject(e);
          });
      });
    },

    errorToMessage(error) {
      if (!error) return "Unknown error";
      if (typeof error === "string") return error;
      if (error.message) return error.message;
      if (error.json && error.json.message) return error.json.message;
      try {
        return JSON.stringify(error);
      } catch (_) {
        return String(error);
      }
    },

    async showRunError(title, result) {
      const detail =
        (result.stderr && result.stderr.trim()) ||
        (result.stdout && result.stdout.trim()) ||
        `Exit code ${result.StatusCode}`;
      await Swal.fire({
        title,
        text: detail,
        confirmButtonText: "OK",
        theme: "dark",
      });
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
