<template>
  <v-tooltip right :disabled="isTooltipDisabled">
    <template v-slot:activator="{ on }">
      <div v-on="on">
        <v-btn
          block
          :disabled="isButtonDisabled"
          :class="buttonClasses"
          @click="handleStartClick"
        >
          Start
        </v-btn>
      </div>
    </template>

    <!-- Tooltip Content -->
    <div v-if="showOptimOTUWarning">
      Missing outgroup database for protax classification
    </div>
    <div v-if="isDockerStopped">
      {{ engineNotFoundMessage }}
    </div>
    <div v-if="isNoFilesSelected">
      No files selected!
    </div>
    <div v-if="showMissingServicesWarning">
      Missing selected services or mandatory inputs
    </div>
    <div v-if="showMissingInputsWarning">
      Missing mandatory inputs
    </div>
  </v-tooltip>
</template>

<script>
import path from 'path';
import os from 'os';
import fs from 'fs';
import slash from 'slash';
import Swal from 'sweetalert2';
import { PassThrough } from 'stream';
import JSONfn from 'json-fn';
import { mapState, mapGetters } from "vuex";
import { stringify } from "envfile";
import cloneDeep from 'lodash/cloneDeep';
import { getServiceScriptsPath } from "../utils/scriptsPath";
import { getContainerUser, prepareBindMounts, applyEngineHostConfig } from "../utils/containerRuntime";

export default {
  name: "Run",
  computed: {
    ...mapState({
      selectedSteps: (state) => state.selectedSteps,
    }),
    ...mapGetters(['isDockerActive', 'engineNotFoundMessage']),
    isButtonDisabled() {
      return this.isDockerStopped ||
             this.isNoFilesSelected ||
             this.isWorkflowNotReady;
    },

    isTooltipDisabled() {
      return !this.isButtonDisabled;
    },

    isDockerStopped() {
      return !this.isDockerActive;
    },

    isNoFilesSelected() {
      return this.$store.state.inputDir === '';
    },

    isWorkflowNotReady() {
      const { workflowName } = this.$route.params;

      if (!workflowName) {
        return !this.$store.getters.selectedStepsReady;
      }

      if (workflowName === 'OptimOTU') {
        return this.$store.state.OptimOTU[8].Inputs[1].value === 'undefined';
      }

      return !this.$store.getters.customWorkflowReady;
    },

    showOptimOTUWarning() {
      return this.$route.params.workflowName === 'OptimOTU' &&
             this.$store.state.OptimOTU[8].Inputs[1].value === 'undefined' ||
             this.$store.state.OptimOTU[8].Inputs[1].value === 'custom';
    },

    showMissingServicesWarning() {
      return !this.$store.getters.selectedStepsReady &&
             !this.$route.params.workflowName;
    },

    showMissingInputsWarning() {
      return 'workflowName' in this.$route.params &&
             !this.$store.getters.customWorkflowReady;
    },

    buttonClasses() {
      return {
        'error-border': this.isButtonDisabled,
        'success-border': !this.isButtonDisabled,
        'bg-dark': true
      };
    }
  },
  data() {
    return {
      userId: null,
      groupId: null
    };
  },
  created() {
    this.initUserAndGroupId();
  },
  methods: {
    initUserAndGroupId() {
      if (os.platform() === 'win32') {
        console.log('Windows system detected, using default user/group IDs');
        this.userId = 1000;
        this.groupId = 1000;
        return;
      }
      try {
        const { execSync } = require('child_process');
        this.userId = parseInt(execSync('id -u').toString().trim());
        this.groupId = parseInt(execSync('id -g').toString().trim());
        console.log(`User ID: ${this.userId}, Group ID: ${this.groupId}`);
      } catch (error) {
        console.warn('Could not get user/group ID, using default');
        this.userId = 1000;
        this.groupId = 1000;
      }
    },
    async confirmRun(name) {
      return Swal.fire({
        title: `Run ${name.replace(/_/g, " ")}?`,
        showCancelButton: true,
        confirmButtonColor: "#3085d6",
        cancelButtonColor: "#d33",
        confirmButtonText: "Continue",
        theme: "dark",
      });
    },
    async updateRunInfo(i, len, Hname, name) {
      this.$store.commit("addRunInfo", [true, name, i, len, Hname]);
    },

    handleStartClick() {
      const { workflowName } = this.$route.params;

      if (!workflowName) {
        return this.runWorkflow({
          name: "workflow",
          steps: this.buildQuickToolSteps(),
        });
      }

      if (workflowName.includes("NextITS")) {
        return this.runWorkflow({
          name: "NextITS",
          steps: this.buildNextITSSteps(),
        });
      }

      if (workflowName.includes("OptimOTU")) {
        return this.runWorkflow({
          name: "OptimOTU",
          steps: this.buildOptimOTUSteps(),
        });
      }

      if (workflowName.includes("FunBarONT")) {
        return this.runWorkflow({
          name: "FunBarONT",
          steps: this.buildFunBarONTSteps(),
          successTitle: "FunBarONT pipeline finished successfully",
          successText: "Results are in your sequences directory",
        });
      }

      return this.runWorkflow({
        name: workflowName,
        steps: this.buildPremadeSteps(workflowName),
      });
    },

    /**
     * One dockerode path for every workflow and every step.
     * Removes the container when it exits. Does not reset workingDir.
     */
    async runContainer(spec) {
      const {
        imageName,
        containerName,
        command,
        env = [],
        binds,
        workingDir,
        log,
        sanitizeChunk,
      } = spec;

      await this.$store.dispatch("imageCheck", imageName);
      await this.$store.dispatch("clearContainerConflicts", containerName);

      const memory = this.$store.state.dockerInfo.MemTotal;
      const nanoCpus = Math.round(Number(this.$store.state.dockerInfo.NCPU) * 1e9);
      const createConfig = {
        Image: imageName,
        name: containerName,
        Cmd: command,
        Tty: false,
        AttachStdout: true,
        AttachStderr: true,
        Platform: "linux/amd64",
        Env: [
          `HOST_UID=${this.userId}`,
          `HOST_GID=${this.groupId}`,
          `fileFormat=${this.$store.state.data.fileFormat}`,
          `readType=${this.$store.state.data.readType}`,
          ...env,
        ],
        HostConfig: applyEngineHostConfig({
          Binds: prepareBindMounts(binds),
          Memory: memory,
          NanoCpus: nanoCpus,
        }),
        User: getContainerUser(this.userId, this.groupId),
      };
      if (workingDir) {
        createConfig.WorkingDir = workingDir;
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
          try { stdoutStream.end(); } catch (err) { console.debug("stdoutStream.end failed:", err && err.message); }
          try { stderrStream.end(); } catch (err) { console.debug("stderrStream.end failed:", err && err.message); }
        };
        attachStream.on("end", endStreams);
        attachStream.on("close", endStreams);

        const logPromise = this.handleDemuxedStreams(
          stdoutStream,
          stderrStream,
          log,
          sanitizeChunk
        );

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
            if (!(msg.includes("HTTP code 404") || msg.includes("HTTP code 409"))) {
              console.warn("Non-fatal remove error:", err);
            }
          }
        }
      }
    },

    /**
     * Shared orchestrator. Premade pipelines pass many steps; OptimOTU,
     * NextITS, FunBarONT, and Quick tools pass one.
     */
    async runWorkflow({ name, steps, successTitle, successText }) {
      let log = null;
      let started = false;
      const startTime = Date.now();

      try {
        const confirmed = await this.confirmRun(name);
        if (!confirmed.isConfirmed) {
          return;
        }
        started = true;

        this.$store.commit("addWorkingDir", "/input");
        this.autoSaveConfig();
        this.$store.state.data.pipeline = name.replace(/ /g, "_");

        if (this.$store.state.data.debugger) {
          log = fs.createWriteStream(
            `${this.$store.state.inputDir}/Pipecraft_${name}_${new Date()
              .toJSON()
              .slice(0, 10)}.txt`
          );
        }

        for (let i = 0; i < steps.length; i++) {
          const step = steps[i];
          if (typeof step.beforeStart === "function") {
            try {
              await step.beforeStart();
            } catch (error) {
              console.error("Failed to generate pipeline configuration:", error);
              await Swal.fire({
                title: "Configuration Error",
                text:
                  error.code === "ENOENT" || error.code === "EACCES"
                    ? "Could not write configuration file. Check file permissions."
                    : error.message || "Failed to generate pipeline configuration.",
                confirmButtonText: "OK",
                theme: "dark",
              });
              return;
            }
          }

          if (step.serviceName) {
            this.$store.state.data.service = String(step.serviceName).replace(
              / /g,
              "_"
            );
          }

          const spec = typeof step.buildSpec === "function"
            ? await step.buildSpec()
            : step;
          const containerName = spec.containerName;
          this.$store.state.runInfo.active = true;
          this.$store.state.runInfo.containerID = containerName;
          this.updateRunInfo(i, steps.length, containerName, name);

          const runResult = await this.runContainer({
            ...spec,
            log,
            sanitizeChunk: spec.sanitizeChunk || step.sanitizeChunk,
          });

          if (runResult.StatusCode === 0) {
            if (typeof step.afterStep === "function") {
              await step.afterStep(runResult);
            }
            this.$store.commit("resetRunInfo");
            continue;
          }

          const message = step.errorFromLogs
            ? this.extractPipelineError(runResult.stdout, runResult.stderr)
            : (runResult.stderr || runResult.stdout || "Unknown error");
          await this.handleDockerError(
            { message, StatusCode: runResult.StatusCode },
            log
          );
          return;
        }

        await Swal.fire({
          title: successTitle || "Workflow finished",
          text: successText || undefined,
          theme: "dark",
        });
      } catch (error) {
        await this.handleDockerError(error, log);
      } finally {
        if (started) {
          this.finishRun(log, startTime);
        }
      }
    },

    finishRun(log, startTime) {
      if (log) {
        try { log.end(); } catch (err) { console.debug("log.end failed:", err && err.message); }
      }
      this.$store.commit("addWorkingDir", "/input");
      this.$store.commit("resetRunInfo");
      if (startTime) {
        console.log(`Total execution time: ${this.toMinsAndSecs(Date.now() - startTime)}`);
      }
    },

    buildPremadeSteps(workflowName) {
      const workflow = this.$store.state[workflowName] || [];
      return workflow
        .filter((step) => step.selected === true || step.selected === "always")
        .map((step) => ({
          serviceName: step.serviceName,
          parseLog: true,
          afterStep: (runResult) => this.applyStepLog(runResult.stdout),
          buildSpec: async () => {
            const dockerProps = await this.getDockerProps(step, workflowName);
            let scriptName = step.scriptName;
            if (typeof scriptName === "object") {
              scriptName = scriptName[this.$store.state.data.dada2mode];
            }
            return {
              imageName: step.imageName,
              containerName: dockerProps.name,
              command: ["bash", "-c", `bash /scripts/${scriptName}`],
              env: dockerProps.Env,
              binds: dockerProps.HostConfig.Binds,
              workingDir: dockerProps.WorkingDir,
            };
          },
        }));
    },

    buildQuickToolSteps() {
      return this.selectedSteps.map((entry, i) => {
        const selectedStep = this.findSelectedService(i);
        return {
          serviceName: selectedStep.serviceName,
          buildSpec: async () => {
            const dockerProps = await this.getDockerProps(selectedStep);
            return {
              imageName: selectedStep.imageName,
              containerName: dockerProps.name,
              command: ["bash", "-c", `bash /scripts/${selectedStep.scriptName}`],
              env: dockerProps.Env,
              binds: dockerProps.HostConfig.Binds,
              workingDir: dockerProps.WorkingDir,
            };
          },
        };
      });
    },

    buildOptimOTUSteps() {
      return [
        {
          serviceName: "optimotu",
          beforeStart: async () => {
            await this.$store.dispatch("generateOptimOTUYamlConfig");
          },
          buildSpec: () => ({
            imageName: "pipecraft/optimotu:5.1-pc1.2.0",
            containerName: "optimotu",
            command: ["/scripts/run_optimotu_dev.sh"],
            env: [
              "R_ENABLE_JIT=0",
              "R_COMPILE_PKGS=0",
              "R_DISABLE_BYTECODE=1",
              "R_KEEP_PKG_SOURCE=yes",
              "LANG=C.UTF-8",
              "LC_ALL=C.UTF-8",
              "LC_CTYPE=C.UTF-8",
              `HOST_OS=${this.$store.state.systemSpecs.os}`,
              `HOST_ARCH=${this.$store.state.systemSpecs.architecture}`,
              `rawFilesDir=${path.basename(this.$store.state.inputDir)}`,
              "R_CLI_NUM_COLORS=0",
              "R_CLI_NO_COLORS=true",
              "NO_COLOR=1",
            ],
            binds: this.getOptimOTUBinds(),
          }),
        },
      ];
    },

    buildFunBarONTSteps() {
      return [
        {
          serviceName: "funbaront",
          errorFromLogs: true,
          beforeStart: async () => {
            await this.$store.dispatch("generateFunBarONTConfig");
          },
          buildSpec: () => ({
            imageName: "pipecraft/funbaront:1-pc1.2.0",
            containerName: "funbaront",
            command: ["/bin/bash", "-c", "bash /scripts/submodules/FunBarONT_Pipeline.sh"],
            env: [
              `HOST_OS=${this.$store.state.systemSpecs.os}`,
              `HOST_ARCH=${this.$store.state.systemSpecs.architecture}`,
              `rawFilesDir=${path.basename(this.$store.state.inputDir)}`,
            ],
            binds: this.getFunBarONTBinds(),
          }),
        },
      ];
    },

    buildNextITSSteps() {
      return [
        {
          serviceName: "Step_1",
          errorFromLogs: true,
          sanitizeChunk: (chunk) => this.sanitizeNextITSLog(chunk),
          beforeStart: async () => {
            await this.$store.dispatch("clearContainerConflicts", "Step_2");
          },
          buildSpec: () => {
            const step = cloneDeep(this.$store.state.NextITS[0]);
            step.Inputs = step.Inputs.concat(this.$store.state.NextITS[1].Inputs);
            step.extraInputs = step.extraInputs.concat(
              this.$store.state.NextITS[1].extraInputs
            );
            const props = this.createParamsFile(step);
            return {
              imageName: "pipecraft/nextits:1.1.0-pc1.2.0",
              containerName: "Step_1",
              command: ["bash", "-c", "bash /scripts/NextITS_Pipeline.sh"],
              env: props.Env,
              binds: props.HostConfig.Binds,
              workingDir: props.WorkingDir,
            };
          },
        },
      ];
    },

    async getDockerProps(step, workflowName) {
      const Hostname = step.serviceName.replaceAll(" ", "_");
      const WorkingDir = this.$store.state.workingDir;
      const envVariables = this.createCustomVariableObj(
        step,
        workflowName || this.$route.params.workflowName
      );
      const Binds = this.getBinds_c(step, this.$store.state.inputDir);
      return {
        Tty: false,
        WorkingDir: WorkingDir,
        User: getContainerUser(this.userId, this.groupId),
        name: Hostname,
        platform: "linux/amd64",
        Volumes: {},
        HostConfig: applyEngineHostConfig({
          Binds: Binds,
          Memory: this.$store.state.dockerInfo.MemTotal,
          NanoCpus: Math.round(Number(this.$store.state.dockerInfo.NCPU) * 1e9)
        }),
        Env: [
          `HOST_UID=${this.userId}`,
          `HOST_GID=${this.groupId}`,
          ...envVariables
        ],
      };
    },

    applyStepLog(stdout) {
      const newWorkingDir = this.getVariableFromLog(stdout, "workingDir");
      const newDataInfo = {
        fileFormat: this.getVariableFromLog(stdout, "fileFormat"),
        output_fasta: this.getVariableFromLog(stdout, "output_fasta"),
        output_feature_table: this.getVariableFromLog(
          stdout,
          "output_feature_table"
        ),
      };
      // Keep data.readType as the workdir choice. Merge/assemble logs
      // print readType=single_end; applying that mid-run would hide the
      // merge step and rewrite PE/SE script names for every workflow.
      this.$store.commit("addInputInfo", {
        fileFormat: newDataInfo.fileFormat || this.$store.state.data.fileFormat,
        readType: this.$store.state.data.readType,
        output_fasta: newDataInfo.output_fasta,
        output_feature_table: newDataInfo.output_feature_table,
      });
      if (newWorkingDir) {
        this.$store.commit("addWorkingDir", newWorkingDir);
      }
    },

    getVariableFromLog(log, varName) {
      try {
        var re = new RegExp(`(${varName}=.*)`, "g");
        const matches = log.match(re);

        if (!matches || matches.length === 0) {
          console.warn(`No match found for variable: ${varName}`);
          return null;
        }

        let value = matches[0].replace('"', "").split("=")[1];
        return value || null;
      } catch (error) {
        console.error(`Error parsing ${varName} from log:`, error);
        return null;
      }
    },
    createCustomVariableObj(element, workflowName) {
      let envVariables = [];
      let nextFlowParams = {};
      let inputs = element.Inputs.concat(element.extraInputs);

      if (Array.isArray(element.extraEnvFromServices)) {
        const resolvedName =
          workflowName || this.$route.params.workflowName;
        const workflow = this.$store.state[resolvedName];
        if (Array.isArray(workflow)) {
          element.extraEnvFromServices.forEach((otherName) => {
            const other = workflow.find((s) => s.serviceName === otherName);
            if (other) {
              inputs = inputs.concat(
                other.Inputs || [],
                other.extraInputs || []
              );
            }
          });
        }
      }

      inputs.forEach((input) => {
        let varObj = {};
        if (input.type === "boolfile") {
          if (input.active === true && input.value != "undefined" && input.value != "") {
            if (Array.isArray(input.value)) {
              nextFlowParams[input.name] = input.value.join();
            } else if (input.name == "ITSx_evalue") {
              nextFlowParams[input.name] = parseFloat(input.value);
            } else if (input.name == "chimera_db") {
              nextFlowParams[input.name] = `/extraFiles15/${path.basename(input.value)}`;
            } else {
              nextFlowParams[input.name] = input.value;
            }
            varObj[input.name] = input.value;
          } else {
            varObj[input.name] = "undefined";
          }
        } else {
          if (input.value != "undefined" && input.value != "") {
            if (Array.isArray(input.value)) {
              nextFlowParams[input.name] = input.value.join();
            } else if (input.name == "ITSx_evalue") {
              nextFlowParams[input.name] = parseFloat(input.value);
            } else if (input.name == "chimera_db") {
              nextFlowParams[input.name] = `/extraFiles15/${path.basename(input.value)}`;
            } else {
              nextFlowParams[input.name] = input.value;
            }
          }
          varObj[input.name] = input.value;
        }
        envVariables.push(stringify(varObj).replace(/(\r\n|\n|\r)/gm, ""));
      });
      let dataInfo = {
        cores: this.$store.state.dockerInfo.NCPU,
        memoryBytes: this.$store.state.dockerInfo.MemTotal,
        workingDir: this.$store.state.workingDir,
        fileFormat: this.$store.state.data.fileFormat,
        readType: this.$store.state.data.readType,
        debugger: this.$store.state.data.debugger,
        dada2mode: this.$store.state.data.dada2mode,
        pipeline: this.$store.state.data.pipeline,
        service: this.$store.state.data.service,
        output_fasta: this.$store.state.data.output_fasta,
        output_feature_table: this.$store.state.data.output_feature_table,
      };
      Object.entries(dataInfo).forEach(([key, value]) => {
        let varObj = {};
        varObj[key] = value;
        envVariables.push(stringify(varObj).replace(/(\r\n|\n|\r)/gm, ""));
      });
      let NextFlowConfigPath = `${getServiceScriptsPath()}/NextFlowConfig.json`;
      if (element.serviceName == "Step_1") {
        fs.writeFile(
          NextFlowConfigPath,
          JSON.stringify(nextFlowParams),
          (error) => {
            if (error) throw error;
          }
        );
      }
      return envVariables;
    },
    getBinds_c(element, Input) {
      const scriptsPath = getServiceScriptsPath();
      let Binds = [`${scriptsPath}:/scripts`, `${Input}:/input`];
      let serviceInputs = element.Inputs.concat(element.extraInputs);
      serviceInputs.forEach((input, index) => {
        if (
          (input.type == "file" &&
            (input.depends_on == undefined ||
              !this.$store.getters.check_depends_on(input))) ||
          (input.type == "boolfile" && input.active == true)
        ) {
          let correctedPath = path.dirname(slash(input.value));
          if (index == 0) {
            let bind = `${correctedPath}:/extraFiles`;
            console.log(bind);
            Binds.push(bind);
          } else {
            let bind = `${correctedPath}:/extraFiles${index + 1}`;
            console.log(bind);
            Binds.push(bind);
          }
        }
      });
      return prepareBindMounts(Binds);
    },
    getOptimOTUBinds() {
      const scriptsPath = getServiceScriptsPath();
      const runsDir = this.$store.state.inputDir;
      let binds = [
        `${scriptsPath}:/scripts`,
        `${runsDir}:/optimotu_targets/sequences`,
        `${runsDir}:/optimotu_targets/sequences/01_raw:rw`,
      ];

      this.$store.state.OptimOTU.forEach((service) => {
        const allInputs = [...(service.Inputs || []), ...(service.extraInputs || [])];
        allInputs.forEach((input) => {
          if (input.type === "boolfile" && input.active === true && input.value) {
            const correctedPath = path.dirname(slash(input.value));

            if (input.name === "custom_sample_table") {
              binds.push(`${correctedPath}:/optimotu_targets/custom_sample_tables`);
            }
            else if (input.name === "positive_control") {
              binds.push(`${correctedPath}:/optimotu_targets/positive_control`);
            }
            else if (input.name === "spike_in") {
              binds.push(`${correctedPath}:/optimotu_targets/spike_in`);
            }
          }
          if (input.name === "cluster_thresholds" &&
              input.value !== "Fungi_GSSP" &&
              input.value !== "Metazoa_MBRAVE") {

            const correctedPath = path.dirname(slash(input.value));
            binds.push(`${correctedPath}:/optimotu_targets/metadata/custom_thresholds`);
          }

          if (input.name === "model_file" &&
              input.value !== "ITS3_ITS4.cm" &&
              input.value !== "f/gITS7_ITS4.cm" &&
              input.value !== "COI.hmm") {

            const correctedPath = path.dirname(slash(input.value));
            binds.push(`${correctedPath}:/optimotu_targets/data/custom_models`);
          }

          if (input.name === "with_outgroup" &&
              input.value !== "UNITE_SHs") {

            const correctedPath = path.dirname(slash(input.value));
            binds.push(`${correctedPath}:/optimotu_targets/data/outgroup`);
          }

          if (input.name === "location" &&
              input.value !== "protaxFungi" &&
              input.value !== "protaxAnimal") {

            const correctedPath = path.dirname(slash(input.value));
            binds.push(`${correctedPath}:/optimotu_targets/protaxCustom`);
          }
        });
      });
      console.log("OptimOTU container binds:", binds);
      return prepareBindMounts(binds);
    },
    getFunBarONTBinds() {
      const taxonomyConfig = this.$store.state.FunBarONT[2];
      const workDir = this.$store.state.inputDir || "";
      const databaseFile = taxonomyConfig?.Inputs?.find(i => i.name === 'database_file')?.value || "";

      if (!databaseFile) {
        throw new Error("No database file selected for FunBarONT (database_file).");
      }

      const scriptDir = getServiceScriptsPath();
      const configPath = `${scriptDir}/FunBarONTConfig.json`;

      return prepareBindMounts([
        `${workDir}:/Input:rw`,
        `${workDir}:/sequences:rw`,
        `${slash(databaseFile)}:/database/database.fasta:ro`,
        `${configPath}:/scripts/FunBarONTConfig.json:ro`,
        `${scriptDir}:/scripts:ro`
      ]);
    },
    findSelectedService(i) {
      let result;
      this.selectedSteps[i].services.forEach((input) => {
        if (input.selected === true || input.selected == "always") {
          result = input;
        }
      });
      return result;
    },
    toMinsAndSecs(millis) {
      var minutes = Math.floor(millis / 60000);
      var seconds = ((millis % 60000) / 1000).toFixed(0);
      return minutes + ":" + (seconds < 10 ? "0" : "") + seconds;
    },
    autoSaveConfig() {
      var conf = [];
      let confJson;
      if (this.$route.params.workflowName) {
        conf.push(this.$store.state[this.$route.params.workflowName]);
        conf.push(this.$route.params.workflowName);
        confJson = JSONfn.stringify(conf);
      } else {
        confJson = JSONfn.stringify(this.$store.state.selectedSteps);
      }
      fs.writeFileSync(
        `${this.$store.state.inputDir}/pipecraft2_last_run_configuration.json`,
        confJson
      );
    },
    createParamsFile(step) {
      let Hostname = step.serviceName.replaceAll(" ", "_");
      let WorkingDir = "/";
      let envVariables = this.createCustomVariableObj(step);
      let Binds = this.getBinds_c(step, this.$store.state.inputDir);
      Binds = prepareBindMounts(Binds.map(b => b.replace(/:\/input(:|$)/, ':/Input$1')));
      return {
        Tty: false,
        WorkingDir: WorkingDir,
        name: Hostname,
        platform: "linux/amd64",
        User: getContainerUser(this.userId, this.groupId),
        Volumes: {},
        HostConfig: applyEngineHostConfig({
          Binds: Binds,
          Memory: this.$store.state.dockerInfo.MemTotal,
          NanoCpus: Math.round(Number(this.$store.state.dockerInfo.NCPU) * 1e9)
        }),
        Env: [
          `HOST_UID=${this.userId}`,
          `HOST_GID=${this.groupId}`,
          ...envVariables
        ],
      };
    },
    sanitizeNextITSLog(chunk) {
      const escChar = String.fromCharCode(27);
      const ansiEscapePattern = new RegExp(
        `${escChar}\\[[0-9;]*[A-Za-z]`,
        "g"
      );
      const stripControlChars = (text) => {
        let result = "";
        for (let i = 0; i < text.length; i += 1) {
          const code = text.charCodeAt(i);
          if (code === 9 || code === 10) {
            result += text[i];
          } else if (code >= 32 && code !== 127) {
            result += text[i];
          }
        }
        return result;
      };
      return stripControlChars(
        chunk
          .replace(ansiEscapePattern, "")
          .replace(/[\u2580-\u259F]/g, "")
          .replace(/\r/g, "")
      );
    },
    handleDemuxedStreams(stdoutStream, stderrStream, log, sanitizeChunk) {
      return new Promise((resolve) => {
        let stdout = '';
        let stderr = '';
        const clean = (text) => (sanitizeChunk ? sanitizeChunk(text) : text);

        const onStdout = (data) => {
          const text = clean(data.toString());
          console.log(text);
          stdout += text;
          if (log) log.write(text);
        };
        const onStderr = (data) => {
          const text = clean(data.toString());
          console.log(text);
          stderr += text;
          if (log) log.write(text);
        };

        stdoutStream.on('data', onStdout);
        stderrStream.on('data', onStderr);

        let endedStdout = false;
        let endedStderr = false;
        const tryResolve = () => {
          if (endedStdout && endedStderr) {
            resolve({ stdout, stderr });
          }
        };

        stdoutStream.on('end', () => { endedStdout = true; tryResolve(); });
        stderrStream.on('end', () => { endedStderr = true; tryResolve(); });
      });
    },

    extractPipelineError(stdout = '', stderr = '') {
      const combined = `${stdout || ''}\n${stderr || ''}`;
      const cleaned = combined
        .split(/\r?\n/)
        .filter((line) => !/Nextflow\s+\S+\s+is available - Please consider updating/i.test(line))
        .join('\n');

      const match = cleaned.match(/(ERROR ~[\s\S]*|[^\n]*input file name collision[\s\S]*|Caused by:[\s\S]*)/);
      const focused = (match ? match[0] : cleaned).trim();
      return focused || 'Unknown error';
    },

    waitWithTimeout(promise, ms) {
      return new Promise((resolve, reject) => {
        const t = setTimeout(() => reject(new Error('timeout')), ms);
        promise.then((v) => { clearTimeout(t); resolve(v); })
               .catch((e) => { clearTimeout(t); reject(e); });
      });
    },

    async handleDockerError(error, log) {
      const statusCode = error?.StatusCode ?? null;
      const message = error?.message || '';

      const isGracefulStop =
        statusCode === 137 ||
        message.includes('HTTP code 404') ||
        message.includes('HTTP code 409');

      const toReadable = (err) => {
        if (!err) return 'Unknown error';
        if (typeof err === 'string') return err;
        if (err.message && typeof err.message === 'string') return err.message;
        try { return JSON.stringify(err); } catch (_) { return String(err); }
      };

      const readable = toReadable(error);

      if (isGracefulStop) {
        console.info('Docker stop detected:', readable);
      } else {
        console.error('Docker error:', readable);
      }
      if (log) {
        log.write(`Error: ${readable}\n`);
      }

      if (isGracefulStop) {
        await Swal.fire({
          title: "Workflow stopped",
          theme: "dark",
        });
        return;
      }

      let extra = '';
      try {
        const dataDir = path.dirname(this.$store.state.inputDir);
        const logPath = path.join(dataDir, 'optimotu_targets.log');
        const content = await fs.promises.readFile(logPath, 'utf8');
        const lines = content.split(/\r?\n/);
        extra = lines.slice(-50).join('\n').trim();
      } catch (_) {
        // ignore if log not available
      }

      const summary = extra && (!readable || readable === 'Unknown error')
        ? extra
        : (extra ? `${readable}\n\n--- Last log lines ---\n${extra}` : readable);

      const esc = (s) => s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');

      await Swal.fire({
        title: "An error has occurred while processing your data",
        html: `<pre style="text-align:left;white-space:pre-wrap;max-height:50vh;overflow:auto">${esc(summary)}</pre>`,
        confirmButtonText: "OK",
        theme: "dark",
        width: 900
      });
    },
  },
};
</script>

<style scoped>
.v-btn {
  justify-content: center;
}

.swal-wide {
  width: 850px !important;
}

.swal2-popup {
  width: auto;
}

.bg-dark {
  background-color: #212121 !important;
}

.error-border {
  border: thin solid #E57373 !important;
  border-top: thin solid white !important;
  border-right: thin solid white !important;
  border-left: thin solid white !important;
}

.success-border {
  border-bottom: thin solid #1DE9B6 !important;
  border-top: thin solid white !important;
  border-left: thin solid white !important;
  border-right: thin solid white !important;
}
</style>
