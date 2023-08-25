<template>
  <v-tooltip
    right
    :disabled="
      $store.state.dockerStatus == 'running' &&
      $store.state.inputDir != '' &&
      (('workflowName' in $route.params &&
        $store.getters.customWorkflowReady) ||
        $store.getters.selectedStepsReady)
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
              !$store.getters.customWorkflowReady)
          "
          block
          :style="
            $store.state.dockerStatus == 'stopped' ||
            $store.state.inputDir == '' ||
            (!('workflowName' in $route.params) &&
              !$store.getters.selectedStepsReady) ||
            ('workflowName' in $route.params &&
              !$store.getters.customWorkflowReady)
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
              : $route.params.workflowName
              ? runCustomWorkFlow($route.params.workflowName)
              : runWorkFlow()
          "
        >
          Start
        </v-btn>
      </div>
    </template>
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
var socketPath =
  os.platform() === "win32" ? "//./pipe/docker_engine" : "/var/run/docker.sock";
var dockerode = new Dockerode({ socketPath: socketPath });
var stdout = new streams.WritableStream();
var stderr = new streams.WritableStream();
const isDevelopment = process.env.NODE_ENV !== "production";
const fs = require("fs");
var JSONfn = require("json-fn");

export default {
  name: "Run",
  computed: mapState({
    selectedSteps: (state) => state.selectedSteps,
  }),
  data: () => ({
    items: [],
  }),
  methods: {
    async confirmRun(name) {
      let result = await Swal.fire({
        title: `Run ${name.replace(/_/g, " ")}?`,
        showCancelButton: true,
        confirmButtonColor: "#3085d6",
        cancelButtonColor: "#d33",
        confirmButtonText: "Continue",
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
        name: Hostname,
        Volumes: {},
        HostConfig: {
          Binds: Binds,
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
              let dockerProps = await this.getDockerProps(step);
              this.updateRunInfo(i, steps2Run, dockerProps.name, name);
              await this.imageCheck(step.imageName);
              await this.clearContainerConflicts(dockerProps.name);
              console.log(dockerProps);
              let scriptName;
              if (typeof step.scriptName === "object") {
                scriptName = step.scriptName[this.$store.state.data.dada2mode];
              } else {
                scriptName = step.scriptName;
              }
              console.log(scriptName);
              let result = await dockerode
                .run(
                  step.imageName,
                  ["sh", "-c", `/scripts/${scriptName}`],
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
                  Swal.fire("Workflow stopped");
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
                  Swal.fire("Workflow finished");
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
          for (let [i, step] of this.selectedSteps.entries()) {
            let selectedStep = this.findSelectedService(i);
            let dockerProps = await this.getDockerProps(selectedStep);
            console.log(dockerProps);
            this.updateRunInfo(i, steps2Run, dockerProps.name, "workflow");
            await this.imageCheck(selectedStep.imageName);
            await this.clearContainerConflicts(dockerProps.name);
            let result = await dockerode
              .run(
                selectedStep.imageName,
                ["sh", "-c", `/scripts/${selectedStep.scriptName}`],
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
                Swal.fire("Workflow stopped");
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
                Swal.fire("Workflow finished");
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
      var re = new RegExp(`(${varName}=.*)`, "g");
      let value = log.match(re)[0].replace('"', "").split("=")[1];
      return value;
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
        nextFlowParams[input.name] = input.value;
        varObj[input.name] = input.value;
        envVariables.push(stringify(varObj).replace(/(\r\n|\n|\r)/gm, ""));
      });
      let dataInfo = {
        workingDir: this.$store.state.workingDir,
        fileFormat: this.$store.state.data.fileFormat,
        readType: this.$store.state.data.readType,
        debugger: this.$store.state.data.debugger,
        dada2mode: this.$store.state.data.dada2mode,
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
      fs.writeFile(
        NextFlowConfigPath,
        JSON.stringify(nextFlowParams),
        (error) => {
          if (error) throw error;
        }
      );
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
          input.type == "file" ||
          (input.type == "boolfile" && input.active == true)
        ) {
          let correctedPath = path.dirname(slash(input.value));
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
    findAndRemoveContainer() {},
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
      fs.writeFileSync(`${this.$store.state.inputDir}/config.json`, confJson);
    },
    createParamsFile(step) {
      let Hostname = step.serviceName.replaceAll(" ", "_");
      let WorkingDir = "/input";
      let envVariables = this.createCustomVariableObj(step);
      let Binds = this.getBinds_c(step, this.$store.state.inputDir);
      let dockerProps = {
        Tty: false,
        WorkingDir: WorkingDir,
        name: Hostname,
        platform: "linux/amd64",
        Volumes: {},
        HostConfig: {
          Binds: Binds,
        },
        Env: envVariables,
      };
      return dockerProps;
    },
    async runNextITS() {
      var writeLog = this.$store.state.data.debugger;
      this.confirmRun("NextITS").then(async (result) => {
        if (result.isConfirmed) {
          this.$store.state.runInfo.active = true;
          let log;
          if (this.$store.state.data.debugger == true) {
            log = fs.createWriteStream("NextITS_log.txt");
          }
          let stdout = new streams.WritableStream();
          // let step = this.$store.state.NextITS[0];
          let step = _.cloneDeep(this.$store.state.NextITS[0]);
          step.Inputs = step.Inputs.concat(this.$store.state.NextITS[1].Inputs);
          step.extraInputs = step.extraInputs.concat(
            this.$store.state.NextITS[1].extraInputs
          );
          let props = this.createParamsFile(step);
          console.log(props);
          await this.clearContainerConflicts("Step_1");
          await this.clearContainerConflicts("Step_2");
          await this.imageCheck("vmikk/nextits:0.5.0");
          let promise = new Promise((resolve, reject) => {
            dockerode
              .run(
                "vmikk/nextits:0.5.0",
                ["sh", "-c", `/scripts/NextITS_Pipeline.sh`],
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
            Swal.fire("Workflow finished");
          } else {
            Swal.fire({
              title: "An error has occured while processing your data",
              text: "unknows error, check result/pipeline_info/execution_report for more info",
              confirmButtonText: "Quit",
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
