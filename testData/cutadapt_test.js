import * as Dockerode from "dockerode";
import os from "os";
var _ = require("lodash");
const streams = require("memory-streams");
const socketPath =
  os.platform() === "win32" ? "//./pipe/docker_engine" : "/var/run/docker.sock";
const dockerode = new Dockerode({ socketPath: socketPath });
const stdout = new streams.WritableStream();
const stderr = new streams.WritableStream();

let dockerProps = {
  Tty: false,
  WorkingDir: "/input",
  name: "cut_primers",
  Volumes: {},
  HostConfig: {
    Binds: [
      "testData/COI_313bp:/input",
      "src/pipecraft-core/service_scripts:/scripts",
    ],
  },
  Env: [
    "forward_primers=GGWACWGGWTGAACWGTWTAYCCYCC",
    "reverse_primers=TANACYTCNGGRTGNCCRAARAAYCA",
    "mismatches=1",
    "min_overlap=21",
    "seqs_to_keep=keep_all",
    "pair_filter=both",
    "cores=1",
    "no_indels=true",
    "workingDir=/input",
    "fileFormat=fastq.gz",
    "readType=paired_end",
    "debugger=false",
    "dada2mode=FORWARD",
  ],
};

async function runWorkFlow() {
  this.$store.commit("addWorkingDir", "/input");
  console.log(dockerProps);
  let result = await dockerode
    .run(
      selectedStep.imageName,
      ["sh", "-c", `/scripts/cut_primers_paired_end_reads.sh`],
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
      return err;
    });
  console.log(result);
  let totalTime = this.toMinsAndSecs(Date.now() - startTime);
  console.log(totalTime);
}

runWorkFlow();
