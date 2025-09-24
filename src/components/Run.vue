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
      Failed to find docker desktop!
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
import { WritableStream } from 'memory-streams';
import { PassThrough } from 'stream';
import JSONfn from 'json-fn';
import { ipcRenderer } from "electron";
import { mapState, mapGetters } from "vuex";
import { stringify } from "envfile";
import cloneDeep from 'lodash/cloneDeep';
var stdout = new WritableStream();
var stderr = new WritableStream();
const isDevelopment = process.env.NODE_ENV !== "production";



export default {
  name: "Run",
  computed: {
    // Keep the mapState for selectedSteps
    ...mapState({
      selectedSteps: (state) => state.selectedSteps,
    }),
    ...mapGetters(['isDockerActive']),
    // Button State
    isButtonDisabled() {
      return this.isDockerStopped || 
             this.isNoFilesSelected || 
             this.isWorkflowNotReady;
    },

    isTooltipDisabled() {
      return !this.isButtonDisabled;
    },

    // Status Checks
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

    // Warning States
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

    // Button Styling
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
        // Windows - use default values
        console.log('Windows system detected, using default user/group IDs');
        this.userId = 1000;
        this.groupId = 1000;
      } else {
        // Linux/macOS
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
      }
    },
    async confirmRun(name) {
      let result = await Swal.fire({
        title: `Run ${name.replace(/_/g, " ")}?`,
        showCancelButton: true,
        confirmButtonColor: "#3085d6",
        cancelButtonColor: "#d33",
        confirmButtonText: "Continue",
        theme: "dark",
      });
      return result;
    },
    async updateRunInfo(i, len, Hname, name) {
      this.$store.commit("addRunInfo", [true, name, i, len, Hname]);
    },
    async getDockerProps(step) {
      let Hostname = step.serviceName.replaceAll(" ", "_");
      let WorkingDir = this.$store.state.workingDir;
      let envVariables = this.createCustomVariableObj(step);
      let Binds = this.getBinds_c(step, this.$store.state.inputDir);
      let dockerProps = {
        Tty: false,
        WorkingDir: WorkingDir,
        User: `${this.userId}:${this.groupId}`,
        name: Hostname,
        Volumes: {},
        HostConfig: {
          Binds: Binds,
          Memory: this.$store.state.dockerInfo.MemTotal,
          NanoCpus: this.$store.state.dockerInfo.NCPU * 1e9,
        },
        Env: envVariables,
      };
      return dockerProps;
    },
    async runCustomWorkFlow(name) {
      this.confirmRun(name).then(async (result) => {
        if (result.isConfirmed) {
          this.$store.commit("addWorkingDir", "/input");
          let startTime = Date.now();
          let steps2Run = this.$store.getters.steps2Run(name);
          this.autoSaveConfig();
          this.$store.state.data.pipeline = name.replace(/ /g, "_");
          let log;
          if (this.$store.state.data.debugger == true) {
            log = fs.createWriteStream(
              `${this.$store.state.inputDir}/Pipecraft_${name}_${new Date()
                .toJSON()
                .slice(0, 10)}.txt`
            );
          }
          for (let [i, step] of this.$store.state[name].entries()) {
            if (step.selected == true || step.selected == "always") {
              this.$store.state.data.service = step.serviceName.replace(
                / /g,
                "_"
              );
              let dockerProps = await this.getDockerProps(step);
              this.updateRunInfo(i, steps2Run, dockerProps.name, name);
              await this.$store.dispatch('imageCheck', step.imageName);
              await this.$store.dispatch('clearContainerConflicts', dockerProps.name);
              console.log(dockerProps);
              let scriptName;
              if (typeof step.scriptName === "object") {
                console.log(step.scriptName[this.$store.state.data.dada2mode]);
                scriptName = step.scriptName[this.$store.state.data.dada2mode];
              } else {
                scriptName = step.scriptName;
              }
              console.log(scriptName);
              let result = await this.$docker
                .run(
                  step.imageName,
                  ["bash", "-c", `bash /scripts/${scriptName}`],
                  [stdout, stderr],
                  dockerProps
                )
                .then(async ([res, container]) => {
                  console.log(stderr.toString());
                  console.log(stdout.toString());
                  res.stdout = stdout.toString();
                  res.stderr = stderr.toString();
                  if (res.StatusCode != 137) {
                    container.remove({ v: true, force: true });
                  }
                  console.log(res);
                  return res;
                })
                .catch((err) => {
                  console.log(err);
                  this.$store.commit("resetRunInfo");
                  return err;
                });
              console.log(result);
              if (result.StatusCode == 0) {
                if (this.$store.state.data.debugger == true) {
                  log.write(result.stdout.toString().replace(/[\n\r]/g, ""));
                }
                let newWorkingDir = this.getVariableFromLog(
                  result.stdout,
                  "workingDir"
                );
                let newDataInfo = {
                  fileFormat: this.getVariableFromLog(
                    result.stdout,
                    "fileFormat"
                  ),
                  readType: this.getVariableFromLog(result.stdout, "readType"),
                  output_fasta: this.getVariableFromLog(
                    result.stdout,
                    "output_fasta"
                  ),
                  output_feature_table: this.getVariableFromLog(
                    result.stdout,
                    "output_feature_table"
                  ),
                };
                this.$store.commit(
                  "toggle_PE_SE_scripts",
                  newDataInfo.readType
                );
                this.$store.commit("addInputInfo", newDataInfo);
                this.$store.commit("addWorkingDir", newWorkingDir);
              } else {
                if (result.StatusCode == 137) {
                  if (this.$store.state.data.debugger == true) {
                    log.write(result.stderr.toString().replace(/[\n\r]/g, ""));
                  }
                  Swal.fire({
                    title: "Workflow stopped",
                    theme: "dark",
                  });
                } else {
                  let err;
                  if (!result.stderr) {
                    if (this.$store.state.data.debugger == true) {
                      log.write(
                        result.stdout.toString().replace(/[\n\r]/g, "")
                      );
                    }
                    err = result;
                  } else {
                    if (this.$store.state.data.debugger == true) {
                      log.write(
                        result.stderr.toString().replace(/[\n\r]/g, "")
                      );
                    }
                    err = result.stderr;
                  }
                  Swal.fire({
                    title: "An error has occured while processing your data",
                    text: err,
                    confirmButtonText: "Quit",
                    theme: "dark",
                  });
                }
                this.$store.commit("resetRunInfo");
                stdout = new WritableStream();
                stderr = new WritableStream();
                break;
              }
              stdout = new WritableStream();
              stderr = new WritableStream();
              console.log(`Finished step ${i + 1}: ${step.serviceName}`);
              this.$store.commit("resetRunInfo");
              if (result.StatusCode == 0) {
                steps2Run -= 1;
                if (steps2Run == 0) {
                  Swal.fire({
                    title: "Workflow finished",
                    theme: "dark",
                  });
                }
              }
            }
          }
          let totalTime = this.toMinsAndSecs(Date.now() - startTime);
          this.$store.commit("addWorkingDir", "/input");
          this.$store.commit("resetRunInfo");
          console.log(totalTime);
        }
      });
    },
    async runWorkFlow() {
      this.confirmRun("workflow").then(async (result) => {
        if (result.isConfirmed) {
          this.$store.commit("addWorkingDir", "/input");
          let startTime = Date.now();
          let steps2Run = this.$store.getters.steps2Run("selectedSteps");
          console.log(`${this.$store.state.inputDir}`);
          this.autoSaveConfig();
          let log;
          if (this.$store.state.data.debugger == true) {
            log = fs.createWriteStream(
              `${
                this.$store.state.inputDir
              }/Pipecraft_CustomWorkflow_${new Date()
                .toJSON()
                .slice(0, 10)}.txt`
            );
          }
          this.$store.state.data.pipeline =
            `quick_tools_${this.selectedSteps.stepName}`.replace(/ /g, "_");
          for (let [i, step] of this.selectedSteps.entries()) {
            let selectedStep = this.findSelectedService(i);

            this.$store.state.data.service = selectedStep.serviceName.replace(
              / /g,
              "_"
            );
            let dockerProps = await this.getDockerProps(selectedStep);
            console.log(dockerProps);
            this.updateRunInfo(i, steps2Run, dockerProps.name, "workflow");
            await this.$store.dispatch('imageCheck', selectedStep.imageName);
            await this.$store.dispatch('clearContainerConflicts', dockerProps.name);
            let result = await this.$docker
              .run(
                selectedStep.imageName,
                ["bash", "-c", `bash /scripts/${selectedStep.scriptName}`],
                [stdout, stderr],
                dockerProps
              )
              .then(async ([res, container]) => {
                res.stdout = stdout.toString();
                res.stderr = stderr.toString();
                if (res.StatusCode != 137) {
                  container.remove({ v: true, force: true });
                }
                console.log(res);
                return res;
              })
              .catch((err) => {
                console.log(err);
                this.$store.commit("resetRunInfo");
                return err;
              });
            console.log(result);
            if (result.StatusCode == 0) {
              if (this.$store.state.data.debugger == true) {
                log.write(result.stdout.toString().replace(/[\n\r]/g, ""));
              }
              let newWorkingDir = this.getVariableFromLog(
                result.stdout,
                "workingDir"
              );
              let newDataInfo = {
                fileFormat: this.getVariableFromLog(
                  result.stdout,
                  "fileFormat"
                ),
                readType: this.getVariableFromLog(result.stdout, "readType"),
              };
              this.$store.commit("addInputInfo", newDataInfo);
              this.$store.commit("addWorkingDir", newWorkingDir);
            } else {
              if (result.StatusCode == 137) {
                if (this.$store.state.data.debugger == true) {
                  log.write(result.stderr.toString().replace(/[\n\r]/g, ""));
                }
                Swal.fire({
                  title: "Workflow stopped",
                  theme: "dark",
                });
              } else {
                let err;
                if (!result.stderr) {
                  if (this.$store.state.data.debugger == true) {
                    log.write(result.stdout.toString().replace(/[\n\r]/g, ""));
                  }
                  err = result;
                } else {
                  err = result.stderr;
                  if (this.$store.state.data.debugger == true) {
                    log.write(result.stderr.toString().replace(/[\n\r]/g, ""));
                  }
                }
                Swal.fire({
                  title: "An error has occured while processing your data",
                  text: err,
                  confirmButtonText: "Quit",
                  theme: "dark",
                });
              }
              this.$store.commit("resetRunInfo");
              stdout = new WritableStream();
              stderr = new WritableStream();
              break;
            }
            stdout = new WritableStream();
            stderr = new WritableStream();
            console.log(`Finished step ${i + 1}: ${step.stepName}`);
            this.$store.commit("resetRunInfo");
            if (result.StatusCode == 0) {
              // steps2Run -= 1;
              if (steps2Run == 0) {
                Swal.fire({
                  title: "Workflow finished",
                  theme: "dark",
                });
              }
            }
          }
          let totalTime = this.toMinsAndSecs(Date.now() - startTime);
          this.$store.commit("addWorkingDir", "/input");
          this.$store.commit("resetRunInfo");
          console.log(totalTime);
        }
      });
    },
    getVariableFromLog(log, varName) {
      try {
        var re = new RegExp(`(${varName}=.*)`, "g");
        const matches = log.match(re);
        
        // Check if we found any matches
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
    createVariableObj(stepIndex, serviceIndex) {
      let envVariables = [];
      this.selectedSteps[stepIndex].services[serviceIndex].Inputs.forEach(
        (input) => {
          let varObj = {};
          varObj[input.name] = input.value;
          envVariables.push(stringify(varObj).replace(/(\r\n|\n|\r)/gm, ""));
        }
      );
      this.selectedSteps[stepIndex].services[serviceIndex].extraInputs.forEach(
        (input) => {
          let varObj = {};
          varObj[input.name] = input.value;
          envVariables.push(stringify(varObj).replace(/(\r\n|\n|\r)/gm, ""));
        }
      );
      let dataInfo = {
        workingDir: this.$store.state.workingDir,
        fileFormat: this.$store.state.data.fileFormat,
        readType: this.$store.state.data.readType,
        debugger: this.$store.sate.data.debugger,
        dada2mode: this.$store.state.data.dada2mode,
        pipeline: this.$store.state.data.pipeline,
        service: this.$store.state.data.service,
      };
      Object.entries(dataInfo).forEach(([key, value]) => {
        let varObj = {};
        varObj[key] = value;
        envVariables.push(stringify(varObj).replace(/(\r\n|\n|\r)/gm, ""));
      });
      return envVariables;
    },
    createCustomVariableObj(element) {
      let envVariables = [];
      let nextFlowParams = {};
      let inputs = element.Inputs.concat(element.extraInputs);
      inputs.forEach((input) => {
        let varObj = {};
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
        envVariables.push(stringify(varObj).replace(/(\r\n|\n|\r)/gm, ""));
      });
      let dataInfo = {
        cores: this.$store.state.dockerInfo.NCPU,
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
      let NextFlowConfigPath =
        isDevelopment == true
          ? `${slash(
              process.cwd()
            )}/src/pipecraft-core/service_scripts/NextFlowConfig.json`
          : `${process.resourcesPath}/src/pipecraft-core/service_scripts/NextFlowConfig.json`;
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
      let scriptsPath =
        isDevelopment == true
          ? `${slash(process.cwd())}/src/pipecraft-core/service_scripts`
          : `${process.resourcesPath}/src/pipecraft-core/service_scripts`;
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
      return Binds;
    },
    createBinds(serviceIndex, stepIndex, Input) {
      let scriptsPath =
        isDevelopment == true
          ? `${slash(process.cwd())}/src/pipecraft-core/service_scripts`
          : `${process.resourcesPath}/src/pipecraft-core/service_scripts`;
      let Binds = [`${scriptsPath}:/scripts`, `${Input}:/input`];
      let serviceInputs = this.selectedSteps[stepIndex].services[
        serviceIndex
      ].Inputs.concat(
        this.selectedSteps[stepIndex].services[serviceIndex].extraInputs
      );
      serviceInputs.forEach((input, index) => {
        if (
          input.type == "file" ||
          (input.type == "boolfile" && input.active == true)
        ) {
          let correctedPath = path.dirname(slash(input.value));
          // let fileName = path.parse(correctedPath).base;
          if (index == 0) {
            let bind = `${correctedPath}:/extraFiles`;
            Binds.push(bind);
          } else {
            let bind = `${correctedPath}:/extraFiles${index + 1}`;
            Binds.push(bind);
          }
        }
      });
      return Binds;
    },
    getOptimOTUBinds() {
      const scriptsPath = isDevelopment 
        ? `${slash(process.cwd())}/src/pipecraft-core/service_scripts`
        : `${process.resourcesPath}/src/pipecraft-core/service_scripts`;
      const runsDir = this.$store.state.inputDir; // user-selected folder containing Run*/
      let binds = [
        `${scriptsPath}:/scripts`,
        // Mount runsDir as sequences root (outputs will be written here)
        `${runsDir}:/optimotu_targets/sequences`,
        // Mount runsDir again as 01_raw (read-only view for inputs)
        `${runsDir}:/optimotu_targets/sequences/01_raw:ro`,
      ];

      // Process all inputs from OptimOTU workflow
      this.$store.state.OptimOTU.forEach((service) => {
        // Combine regular and extra inputs
        const allInputs = [...(service.Inputs || []), ...(service.extraInputs || [])];
        allInputs.forEach((input) => {
          // Handle specific boolfile inputs
          if (input.type === "boolfile" && input.active === true && input.value) {
            const correctedPath = path.dirname(slash(input.value));

            // Handle each boolfile input specifically
            if (input.name === "custom_sample_table") {
              binds.push(`${correctedPath}:/optimotu_targets/custom_sample_tables`);
              // Note: In your YAML, you'd reference this as /optimotu_targets/data/sample_tables/${fileName}
            } 
            else if (input.name === "positive_control") {
              binds.push(`${correctedPath}:/optimotu_targets/positive_control`);
              // Note: In your YAML, you'd reference this as /optimotu_targets/data/controls/positive/${fileName}
            }
            else if (input.name === "spike_in") {
              binds.push(`${correctedPath}:/optimotu_targets/spike_in`);
              // Note: In your YAML, you'd reference this as /optimotu_targets/data/controls/spike_in/${fileName}
            }
          }
          // Special handling for specific inputs
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
      return binds;
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
    async runStep(envVariables, scriptName, imageName) {
      var result = await ipcRenderer.sendSync(
        "runStep",
        imageName,
        scriptName,
        envVariables,
        this.$store.state.workingDir
      );
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
      let WorkingDir = "/input";
      let envVariables = this.createCustomVariableObj(step);
      let Binds = this.getBinds_c(step, this.$store.state.inputDir);
      
      // Replace the /input bind with /input/Input
      const inputBindIndex = Binds.findIndex(bind => bind.includes(':/input'));
      if (inputBindIndex !== -1) {
        Binds[inputBindIndex] = Binds[inputBindIndex].replace(':/input', ':/input/Input');
      }
      
      let dockerProps = {
        Tty: false,
        WorkingDir: WorkingDir,
        User: `${this.userId}:${this.groupId}`,
        name: Hostname,
        platform: "linux/amd64",
        Volumes: {},
        HostConfig: {
          Binds: Binds,
          Memory: this.$store.state.dockerInfo.MemTotal,
          NanoCpus: this.$store.state.dockerInfo.NCPU * 1e9,
        },
        Env: envVariables,
      };
      return dockerProps;
    },
    async runOptimOTU() {
      let container = null;
      let log = null;
      let startTime = null;
      
      try {
        const result = await this.confirmRun('OptimOTU');
        if (!result.isConfirmed) return;
      
        const setup = await this.setupWorkflow('OptimOTU');
        startTime = setup.startTime;
        log = setup.log;
      
        this.$store.state.runInfo.active = true;
        this.$store.state.runInfo.containerID = 'optimotu';
      
        try {
          await this.$store.dispatch('generateOptimOTUYamlConfig');
        } catch (error) {
          console.error('Failed to generate YAML config:', error);
          
          await Swal.fire({
            title: "Configuration Error",
            text: error.code === 'ENOENT' || error.code === 'EACCES' 
              ? "Could not write configuration file. Check file permissions."
              : "Failed to generate pipeline configuration.",
            confirmButtonText: "OK",
            theme: "dark",
          });
          
          this.$store.commit("resetRunInfo");
          return;
        }
      
        const { container: dockerContainer, stdoutStream, stderrStream } = await this.executeDockerContainer({
          imageName: 'pipecraft/optimotu:5',
          containerName: 'optimotu',
          command: ['/scripts/run_optimotu.sh'],
          env: [
            'R_ENABLE_JIT=0',
            'R_COMPILE_PKGS=0',
            'R_DISABLE_BYTECODE=1',
            'R_KEEP_PKG_SOURCE=yes',
            'LANG=C.UTF-8',
            'LC_ALL=C.UTF-8',
            'LC_CTYPE=C.UTF-8',
            `HOST_OS=${this.$store.state.systemSpecs.os}`,
            `HOST_ARCH=${this.$store.state.systemSpecs.architecture}`,
            `rawFilesDir=${path.basename(this.$store.state.inputDir)}`,
            'R_CLI_NUM_COLORS=0',
            'R_CLI_NO_COLORS=true',
            'NO_COLOR=1'
          ],
          binds: this.getOptimOTUBinds(),
          memory: this.$store.state.dockerInfo.MemTotal,
          cpuCount: this.$store.state.dockerInfo.NCPU,
          userId: this.userId,
          groupId: this.groupId
        });
      
        container = dockerContainer;

        // Start log capture but don't block on it yet
        const logPromise = this.handleDemuxedStreams(stdoutStream, stderrStream, log);

        // Wait for the container process to finish (authoritative)
        const data = await container.wait();

        // Ensure streams are closed so logPromise can resolve
        try { stdoutStream.end(); } catch (err) { console.debug('stdoutStream.end failed (likely closed):', err && err.message ? err.message : err); }
        try { stderrStream.end(); } catch (err) { console.debug('stderrStream.end failed (likely closed):', err && err.message ? err.message : err); }

        // Drain any remaining logs with a timeout safeguard
        let stdout = '';
        let stderr = '';
        try {
          const res = await this.waitWithTimeout(logPromise, 2000);
          stdout = res.stdout || '';
          stderr = res.stderr || '';
        } catch (err) {
          console.debug('log drain timeout or error:', err && err.message ? err.message : err);
        }
        await this.cleanupWorkflow(container, log, startTime);
        if (data.StatusCode === 0) {
          await Swal.fire({
            title: "Workflow finished",
            theme: "dark",
          });
        } else {
          await this.handleDockerError({ 
            message: stderr || stdout || 'Unknown error',
            StatusCode: data.StatusCode 
          }, log, container, startTime);
        }
      
      } catch (error) {
        await this.handleDockerError(error, log, container, startTime);
      }
    },
    async generateAndSaveYaml() {
      try {
        const yamlString = await this.$store.dispatch('generateOptimOTUYamlConfig');
        console.log('Generated YAML:', yamlString);
        console.log('YAML configuration generated successfully');
        // Maybe show a success message to the user
      } catch (error) {
        console.error('Error generating YAML configuration:', error);
        // Handle the error, maybe show an error message to the user
      }
    },
    async runNextITS() {
      this.autoSaveConfig();
      var writeLog = this.$store.state.data.debugger;
      this.confirmRun("NextITS").then(async (result) => {
        if (result.isConfirmed) {
          this.$store.state.runInfo.active = true;
          this.$store.state.runInfo.containerID = "Step_1";
          let log;
          if (this.$store.state.data.debugger == true) {
            log = fs.createWriteStream("NextITS_log.txt");
          }
          let stdout = new WritableStream();
          let step = cloneDeep(this.$store.state.NextITS[0]);
          step.Inputs = step.Inputs.concat(this.$store.state.NextITS[1].Inputs);
          step.extraInputs = step.extraInputs.concat(
            this.$store.state.NextITS[1].extraInputs
          );
          let props = this.createParamsFile(step);
          console.log(props);
          await this.$store.dispatch('clearContainerConflicts', "Step_1");
          await this.$store.dispatch('clearContainerConflicts', "Step_2");
          await this.$store.dispatch('imageCheck', "pipecraft/nextits:1.0.0");
          let promise = new Promise((resolve, reject) => {
            this.$docker
              .run(
                "pipecraft/nextits:1.0.0",
                ["bash", "-c", `bash /scripts/NextITS_Pipeline.sh`],
                false,
                props,
                (err, data, container) => {
                  console.log(container);
                  console.log(data);
                  console.log(stdout.toString());
                  if (err) {
                    console.log(err);
                    reject(err);
                  } else {
                    resolve(data);
                  }
                }
              )

              .on("stream", (stream) => {
                stream.on("data", function (data) {
                  console.log(data.toString().replace(/[\n\r]/g, ""));
                  if (writeLog == true) {
                    log.write(data.toString().replace(/[\n\r]/g, ""));
                  }
                  // term.write(data.toString().replace(/[\n\r]/g, "") + "\n");
                  stdout.write(data.toString().replace(/[\n\r]/g, "") + "\n");
                });
              });
          });
          let result = await promise;
          console.log(result);
          this.$store.commit("resetRunInfo");
          if (result.StatusCode == 0) {
            Swal.fire({
              title: "Workflow finished",
              theme: "dark",
            });
          } else {
            Swal.fire({
              title: "An error has occured while processing your data",
              text: "unknown error, check result/pipeline_info/execution_report for more info",
              confirmButtonText: "Quit",
              theme: "dark",
            });
          }
        }
      });
    },
    handleStartClick() {
      const { workflowName } = this.$route.params;
      
      if (!workflowName) {
        return this.runWorkFlow();
      }

      if (workflowName.includes('NextITS')) {
        return this.runNextITS();
      }

      if (workflowName.includes('OptimOTU')) {
        return this.runOptimOTU();
      }

      return this.runCustomWorkFlow(workflowName);
    },
    // Common setup for all workflows
    async setupWorkflow(name) {
      this.$store.commit("addWorkingDir", "/input");
      const startTime = Date.now();
      this.autoSaveConfig();
      
      let log = null;
      if (this.$store.state.data.debugger) {
        log = fs.createWriteStream(
          `${this.$store.state.inputDir}/Pipecraft_${name}_${new Date().toJSON().slice(0, 10)}.txt`
        );
      }
      
      return { startTime, log };
    },

    // Common Docker execution
    async executeDockerContainer(config) {
      const {
        imageName,
        containerName,
        command,
        env,
        binds,
        memory,
        cpuCount,
        userId,
        groupId
      } = config;

      // Check image and clear conflicts
      await this.$store.dispatch('imageCheck', imageName);
      await this.$store.dispatch('clearContainerConflicts', containerName);

      // Create container
      const container = await this.$docker.createContainer({
        Image: imageName,
        name: containerName,
        Cmd: command,
        Tty: false,
        AttachStdout: true,
        AttachStderr: true,
        Platform: "linux/amd64",
        Env: [
          `HOST_UID=${userId}`,
          `HOST_GID=${groupId}`,
          `fileFormat=${this.$store.state.data.fileFormat}`,
          `readType=${this.$store.state.data.readType}`,
          ...env
        ],
        HostConfig: {
          Binds: binds,
          Memory: memory,
          NanoCpus: cpuCount * 1e9,
        }
      });

      // Attach to container (multiplexed stream)
      const attachStream = await container.attach({
        stream: true,
        stdout: true,
        stderr: true,
      });

      // Demux stdout/stderr into separate streams
      const stdoutStream = new PassThrough();
      const stderrStream = new PassThrough();
      // dockerode exposes modem for demuxing multiplexed streams
      // In TTY=false mode, attach() returns multiplexed stream
      container.modem.demuxStream(attachStream, stdoutStream, stderrStream);

      // Propagate end/close so readers resolve
      const endStreams = () => {
        try { stdoutStream.end(); } catch (err) { console.debug('stdoutStream.end on attach close failed:', err && err.message ? err.message : err); }
        try { stderrStream.end(); } catch (err) { console.debug('stderrStream.end on attach close failed:', err && err.message ? err.message : err); }
      };
      attachStream.on('end', endStreams);
      attachStream.on('close', endStreams);

      // Start container
      await container.start();
      
      return { container, stdoutStream, stderrStream };
    },

    // Properly handle demuxed stdout/stderr streams
    handleDemuxedStreams(stdoutStream, stderrStream, log) {
      return new Promise((resolve) => {
        let stdout = '';
        let stderr = '';

        const onStdout = (data) => {
          const text = data.toString();
          console.log(text);
          stdout += text;
          if (log) log.write(text);
        };
        const onStderr = (data) => {
          const text = data.toString();
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

    // Utility: await a promise with timeout
    waitWithTimeout(promise, ms) {
      return new Promise((resolve, reject) => {
        const t = setTimeout(() => reject(new Error('timeout')), ms);
        promise.then((v) => { clearTimeout(t); resolve(v); })
               .catch((e) => { clearTimeout(t); reject(e); });
      });
    },

    // Common error handling
    async handleDockerError(error, log, container = null, startTime = null) {
      const statusCode = error?.StatusCode ?? null;
      const message = error?.message || '';

      // Prefer info logs for expected stop scenarios
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

      // Log
      if (isGracefulStop) {
        console.info('Docker stop detected:', readable);
      } else {
        console.error('Docker error:', readable);
      }
      if (log) {
        log.write(`Error: ${readable}\n`);
      }

      // UX
      if (isGracefulStop) {
        await Swal.fire({
          title: "Workflow stopped",
          theme: "dark",
        });
      } else {
        // Try to append tail of pipeline log for more context
        let extra = '';
        try {
          const dataDir = path.dirname(this.$store.state.inputDir);
          const logPath = path.join(dataDir, 'optimotu_targets.log');
          const content = await fs.promises.readFile(logPath, 'utf8');
          const lines = content.split(/\r?\n/);
          const tail = lines.slice(-50).join('\n');
          extra = tail.trim();
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
      }

      // Always attempt cleanup even on error
      if (container || startTime) {
        await this.cleanupWorkflow(container, log, startTime);
      }
    },

    // Common cleanup (idempotent)
    async cleanupWorkflow(container, log, startTime) {
      if (container) {
        // Try to stop first (ignore if already stopped/removed)
        try {
          await container.stop({ t: 5 });
        } catch (err) {
          // Ignore common non-fatal stop errors
          const msg = err?.message || '';
          if (!(msg.includes('HTTP code 304') || msg.includes('HTTP code 404') || msg.includes('is already in progress'))) {
            console.warn('Non-fatal stop error:', err);
          }
        }

        // Then remove (ignore if already removed or removing)
        try {
          await container.remove({ v: true, force: true });
        } catch (err) {
          const msg = err?.message || '';
          if (!(msg.includes('HTTP code 404') || msg.includes('HTTP code 409'))) {
            console.warn('Non-fatal remove error:', err);
          }
        }
      }
      if (log) {
        log.end();
      }
      this.$store.commit("addWorkingDir", "/input");
      this.$store.commit("resetRunInfo");
      
      if (startTime) {
        const totalTime = this.toMinsAndSecs(Date.now() - startTime);
        console.log(`Total execution time: ${totalTime}`);
      }
    }
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
