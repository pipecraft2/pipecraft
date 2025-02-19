<template>
  <v-card
    light
    elevation="2"
    :disabled="
      Object.values(inputData).includes(input.disabled) ||
      $store.state.runInfo.active == true ||
      $store.getters.check_depends_on(input)
    "
  >
    <v-tooltip top>
      <template v-slot:activator="{ on }">
        <v-card-title
          v-if="$store.getters.linkify(input.tooltip) !== null"
          v-on="on"
          style="justify-content: center; padding: 10px 0px"
        >
          <a :href="$store.getters.linkify(input.tooltip)" target="_blank">{{
            input.name.replace(/_/g, " ")
          }}</a></v-card-title
        >
        <v-card-title
          v-else
          v-on="on"
          style="justify-content: center; padding: 10px 0px"
          >{{ input.name.replace(/_/g, " ") }}</v-card-title
        >
      </template>
      <span>{{ input.tooltip }}</span>
    </v-tooltip>
    <v-card-actions style="justify-content: center">
      <v-row style="justify-content: center"
        ><v-col style="padding: 0" cols="8" offset="0">
          <v-select
            v-model="input.value"
            style="padding-top: 10%"
            @change="inputUpdate(input.value)"
            :items="input.items"
            outlined
          ></v-select>
        </v-col>
      </v-row>
    </v-card-actions>
  </v-card>
</template>

<script>
const slash = require("slash");
const { dialog } = require("@electron/remote");
export default {
  computed: {
    input() {
      if (this.$route.params.workflowName) {
        return this.$store.state[this.$route.params.workflowName][
          this.$attrs.serviceIndex
        ][this.$attrs.list][this.$attrs.inputIndex];
      } else {
        return this.$store.state.selectedSteps[this.$route.params.order]
          .services[this.$attrs.serviceIndex][this.$attrs.list][
          this.$attrs.inputIndex
        ];
      }
    },
    inputData() {
      return this.$store.state.data;
    },
  },
  methods: {
    inputUpdate(value) {
      if (this.$route.params.workflowName) {
        if (value =='custom') {
          this.fileSelect();
        }
        this.blastSwitch(value);
        this.unoiseSwitch(this.$route.params.workflowName, value);
        this.$store.commit("premadeInputUpdate", {
          workflowName: this.$route.params.workflowName,
          serviceIndex: this.$attrs.serviceIndex,
          inputIndex: this.$attrs.inputIndex,
          listName: this.$attrs.list,
          value: value,
        });
      } else {
        if (value =='custom') {
          this.fileSelect();
        }
        this.blastSwitch2(
          value,
          this.$route.params.order,
          this.$attrs.serviceIndex
        );
        this.$store.commit("inputUpdate", {
          stepIndex: this.$route.params.order,
          serviceIndex: this.$attrs.serviceIndex,
          inputIndex: this.$attrs.inputIndex,
          listName: this.$attrs.list,
          value: value,
        });
      }
    },
    blastSwitch(value) {
      if (value == "blastn" || value == "megablast") {
        this.$store.commit("blastSwitch", value);
      }
    },
    blastSwitch2(value, i1, i2) {
      if (value == "blastn" || value == "megablast") {
        this.$store.commit("blastSwitch2", { value: value, i1: i1, i2: i2 });
      }
    },
    unoiseSwitch(name, value) {
      if (name == "NextITS" && value == "unoise") {
        this.$store.state.NextITS[1].Inputs[4].value = true;
      }
      if (name == "NextITS" && (value == "vsearch" || value == "swarm")) {
        this.$store.state.NextITS[1].Inputs[4].value = false;
      }
    },
    fileSelect() {
      dialog
        .showOpenDialog({
          title: "Select input files",
          properties: ["openFile", "multiSelections", "showHiddenFiles"],
        })
        .then((result) => {
          console.log(result);
          if (typeof result.filePaths[0] !== "undefined") {
            let correctedPath = slash(result.filePaths[0]);
            this.inputUpdate(correctedPath);
          }
        })
        .catch((err) => {
          console.log(err);
        });
    },
  },
};
</script>

<style lang="scss" scoped>
.v-text-field {
  ::v-deep input {
    text-align: center !important;
  }
}
.v-text-field input {
  text-align: center;
}
</style>
