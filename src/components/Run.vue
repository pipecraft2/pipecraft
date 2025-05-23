<template>
  <v-tooltip
    right
    :disabled="
      $store.state.dockerStatus == 'running' &&
      $store.state.inputDir != '' &&
      (('workflowName' in $route.params &&
        $store.getters.customWorkflowReady) ||
        $store.getters.selectedStepsReady) &&
    !($route.params.workflowName === 'OptimOTU' && 
      $store.state.OptimOTU[8].Inputs[1].value === 'undefined')
    "
  >
    <template v-slot:activator="{ on }">
      <div v-on="on">
        <v-btn
          style="background-color: #212121"
          :disabled="
            $store.state.dockerStatus == 'stopped' ||
            $store.state.inputDir == '' ||
            (!('workflowName' in $route.params) &&
              !$store.getters.selectedStepsReady) ||
            ('workflowName' in $route.params &&
              !$store.getters.customWorkflowReady) ||
            ($route.params.workflowName === 'OptimOTU' && 
              $store.state.OptimOTU[8].Inputs[1].value === 'undefined')
          "
          block
          :style="
            $store.state.dockerStatus == 'stopped' ||
            $store.state.inputDir == '' ||
            (!('workflowName' in $route.params) &&
              !$store.getters.selectedStepsReady) ||
            ('workflowName' in $route.params &&
              !$store.getters.customWorkflowReady) ||
            ($route.params.workflowName === 'OptimOTU' && 
              $store.state.OptimOTU[8].Inputs[1].value === 'undefined')
              ? {
                  borderBottom: 'thin #E57373 solid',
                  borderTop: 'thin white solid',
                  borderRight: 'thin white solid',
                  borderLeft: 'thin white solid',
                }
              : {
                  borderBottom: 'thin #1DE9B6 solid',
                  borderTop: 'thin white solid',
                  borderLeft: 'thin white solid',
                  borderRight: 'thin white solid',
                }
          "
          @click="
            $route.params.workflowName &&
            $route.params.workflowName.includes('NextITS')
              ? runNextITS()
              : $route.params.workflowName &&
                $route.params.workflowName.includes('OptimOTU')
              ? runOptimOTU()
              : $route.params.workflowName
              ? runCustomWorkFlow($route.params.workflowName)
              : runWorkFlow()
          "
        >
          Start
        </v-btn>
      </div>
    </template>
    <div
    v-if="$route.params.workflowName === 'OptimOTU' && 
          $store.state.OptimOTU[8].Inputs[1].value === 'undefined'"
  >
    Missing outgroup database for protax classification
  </div>
    <div v-if="this.$store.state.dockerStatus == 'stopped'">
      Failed to find docker desktop!
    </div>
    <div v-if="this.$store.state.inputDir == ''">No files selected!</div>
    <div
      v-if="
        !$store.getters.selectedStepsReady &&
        $route.params.workflowName == undefined
      "
    >
      Missing selected services or mandatory inputs
    </div>
    <div
      v-if="
        'workflowName' in $route.params && !$store.getters.customWorkflowReady
      "
    >
      Missing mandatory inputs
    </div>
  </v-tooltip>
</template>

<script>
const path = require("path");
const slash = require("slash");
const Swal = require("sweetalert2");
const streams = require("memory-streams");
var _ = require("lodash");
import * as Dockerode from "dockerode";
import { pullImageAsync } from "dockerode-utils";
import { imageExists } from "dockerode-utils";
import { ipcRenderer } from "electron";
import { mapState } from "vuex";
import { stringify } from "envfile";
import os from "os";
const fs = require('fs');
var JSONfn = require("json-fn");
var socketPath =
  os.platform() === "win32" ? "//./pipe/docker_engine" : "/var/run/docker.sock";
var dockerode = new Dockerode({ socketPath: socketPath });
var stdout = new streams.WritableStream();
var stderr = new streams.WritableStream();
const isDevelopment = process.env.NODE_ENV !== "production";


export default {
  name: "Run",
  computed: mapState({
    selectedSteps: (state) => state.selectedSteps,
  }),
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
    async clearContainerConflicts(Hostname) {
      let container = await dockerode.getContainer(Hostname);
      let nameConflicts = await container
        .remove({ force: true })
        .then(async () => {
          return "Removed conflicting duplicate container";
        })
        .catch(() => {
          return "No conflicting container names";
        });
      console.log(nameConflicts);
    },
    async updateRunInfo(i, len, Hname, name) {
      this.$store.commit("addRunInfo", [true, name, i, len, Hname]);
    },
    async imageCheck(imageName) {
      let gotImg = await imageExists(dockerode, imageName);
      if (gotImg === false) {
        this.$store.commit("activatePullLoader");
        console.log(`Pulling image ${imageName}`);
        let output = await pullImageAsync(dockerode, imageName);
        console.log(output);
        console.log(`Pull complete`);
        this.$store.commit("deactivatePullLoader");
      }
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
              await this.imageCheck(step.imageName);
              await this.clearContainerConflicts(dockerProps.name);
              console.log(dockerProps);
              let scriptName;
              if (typeof step.scriptName === "object") {
                console.log(step.scriptName[this.$store.state.data.dada2mode]);
                scriptName = step.scriptName[this.$store.state.data.dada2mode];
              } else {
                scriptName = step.scriptName;
              }
              console.log(scriptName);
              let result = await dockerode
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
                stdout = new streams.WritableStream();
                stderr = new streams.WritableStream();
                break;
              }
              stdout = new streams.WritableStream();
              stderr = new streams.WritableStream();
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
            await this.imageCheck(selectedStep.imageName);
            await this.clearContainerConflicts(dockerProps.name);
            let result = await dockerode
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
              stdout = new streams.WritableStream();
              stderr = new streams.WritableStream();
              break;
            }
            stdout = new streams.WritableStream();
            stderr = new streams.WritableStream();
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
      const DataDir = this.$store.state.inputDir;
      let binds = [
        `${scriptsPath}:/scripts`,
        `${DataDir}:/optimotu_targets/sequences`,
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
      this.confirmRun('OptimOTU').then(async (result) => {
        if (result.isConfirmed) {
          console.log(this.$route.params.workflowName);
          this.autoSaveConfig();
          this.$store.state.runInfo.active = true;
          this.$store.state.runInfo.containerID = 'optimotu';
          await this.$store.dispatch('generateOptimOTUYamlConfig');
          await this.imageCheck('pipecraft/optimotu:5');
          await this.clearContainerConflicts('optimotu');
          try {            
            const container = await dockerode.createContainer({
              Image: 'pipecraft/optimotu:5',
              name: 'optimotu',
              Cmd: ['/scripts/run_optimotu.sh'],
              Tty: false,
              AttachStdout: true,
              AttachStderr: true,
              Platform: "linux/amd64",
              Env: [
                `HOST_UID=${this.userId}`,
                `HOST_GID=${this.groupId}`,
                `fileFormat=${this.$store.state.data.fileFormat}`,
                `readType=${this.$store.state.data.readType}`,
                'R_ENABLE_JIT=0',                    // Disable JIT compilation
                'R_COMPILE_PKGS=0',                  // Disable package compilation
                'R_DISABLE_BYTECODE=1',              // Disable bytecode compilation
                'R_KEEP_PKG_SOURCE=yes'              // Keep package sources
              ],
              HostConfig: {
                Binds: this.getOptimOTUBinds(),
                Memory: this.$store.state.dockerInfo.MemTotal,
                NanoCpus: this.$store.state.dockerInfo.NCPU * 1e9,
              }
            });
    
            const stream = await container.attach({
              stream: true,
              stdout: true,
              stderr: true,
            });
    
            stream.on('data', (data) => {
              console.log(data.toString());
            });
            
            await container.start();
    
            const data = await container.wait();
            console.log('Container exited with status code:', data.StatusCode);
            this.$store.commit("resetRunInfo");

            if (data.StatusCode == 0) {
              Swal.fire({
                title: "Workflow finished",
                theme: "dark",
              });
            } else {
              Swal.fire({
                title: "An error has occured while processing your data",
                text: "Check the log file for more information",
                confirmButtonText: "Quit",
                theme: "dark",
              });
            }
            await container.remove({ v: true, force: true });
          } catch (err) {
            console.error('Error running container:', err);
            this.$store.commit("resetRunInfo");
            
            // Check if the error is due to the user stopping the container
            if (err.message && err.message.includes('HTTP code 404')) {
              Swal.fire({
                title: "Workflow stopped",
                confirmButtonText: "Quit",
                theme: "dark",
              });
            } else {
              Swal.fire({
                title: "An error has occurred while processing your data",
                text: err.toString(),
                confirmButtonText: "Quit",  
                theme: "dark",
              });
            }
          } 
        }
      });
    },
    generateAndSaveYaml() {
      try {
        const yamlString = this.$store.dispatch('generateOptimOTUYamlConfig');
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
          let stdout = new streams.WritableStream();
          let step = _.cloneDeep(this.$store.state.NextITS[0]);
          step.Inputs = step.Inputs.concat(this.$store.state.NextITS[1].Inputs);
          step.extraInputs = step.extraInputs.concat(
            this.$store.state.NextITS[1].extraInputs
          );
          let props = this.createParamsFile(step);
          console.log(props);
          await this.clearContainerConflicts("Step_1");
          await this.clearContainerConflicts("Step_2");
          await this.imageCheck("pipecraft/nextits:1.0.0");
          let promise = new Promise((resolve, reject) => {
            dockerode
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
</style>
