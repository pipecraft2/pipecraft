<template>
  <!-- each sub component must recive -->
  <v-row justify="center" style="padding-top: 20px">
    <v-card
      height="fit-content"
      width="100%"
      class="mr-4 ml-4 mb-2"
      dark
      elevation="2"
    >
      <v-tooltip top>
        <template v-slot:activator="{ on }">
          <v-card-title
            v-on="on"
            style="justify-content: center; padding: 10px 0px"
            >{{
              $store.state.customWorkflowInfo[$route.params.workflowName].title
            }}</v-card-title
          >
        </template>
        <span></span>
      </v-tooltip>
      <v-card-text
        class="pr-15 pl-15 text-center"
        style="justify-content: center"
        >{{
          this.$store.state.customWorkflowInfo[$route.params.workflowName].info
        }}<br />
        <a
          :href="
            this.$store.state.customWorkflowInfo[$route.params.workflowName]
              .link
          "
          target="_blank"
          >{{
            this.$store.state.customWorkflowInfo[$route.params.workflowName]
              .link
          }}</a
        ></v-card-text
      >
      <v-card-actions v-if="$route.params.workflowName.includes('DADA2')">
        <v-tabs :value="$store.getters.dada2modeIndex">
          <v-tab
            :disabled="
              $store.state.data.readType == 'single_end' &&
              $store.state.inputDir != ''
            "
            @click="$store.commit('setDADAmode', 'FORWARD')"
            >PAIRED-END FORWARD</v-tab
          >
          <v-tab
            :disabled="
              $store.state.data.readType == 'single_end' &&
              $store.state.inputDir != ''
            "
            @click="$store.commit('setDADAmode', 'MIXED')"
            >PAIRED-END MIXED</v-tab
          >
          <v-tab
            :disabled="$store.state.data.readType == 'paired_end'"
            @click="$store.commit('setDADAmode', 'SINGLE_END')"
            >SINGLE-READ</v-tab
          >
        </v-tabs>
      </v-card-actions>
    </v-card>
    <v-expansion-panels
      dark
      multiple
      popout
      style
      :value="
        this.$store.state[this.$route.params.workflowName].length == 1
          ? panel_single
          : panel
      "
    >
      <v-expansion-panel
        v-for="(service, index) in services"
        :key="index"
        :disabled="
          Object.values(inputData).includes(service.disabled) &&
          $store.state.runInfo.active == false
        "
        :class="Object.values(inputData).includes(service.disabled) && hide"
      >
        <v-tooltip top>
          <template v-slot:activator="{ on }">
            <v-expansion-panel-header
              v-on="on"
              style="justify-content: left"
              :class="[service.selected]"
              :disabled="Object.values(inputData).includes(service.disabled)"
            >
              <v-checkbox
                :disabled="Object.values(inputData).includes(service.disabled) || $store.state.runInfo.active == true"
                v-if="service.selected != 'always'"
                hide-details="true"
                @change="check_one($event, index)"
                @click.stop
                v-model="service.selected"
                style="max-width: 34px; padding-top: 0; margin: 0"
              ></v-checkbox>
              {{ service.serviceName.toUpperCase() }}

              <div
                v-if="
                  $store.state.runInfo.active == true &&
                  $store.state.runInfo.type == $route.params.workflowName &&
                  index == $store.state.runInfo.step
                "
              >
                <span
                  v-if="$store.state.pullLoader.active == true"
                  style="padding-left: 25px"
                  >Pulling image ...</span
                >
              </div>
              <div v-if="service.manualLink">
                <v-tooltip left>
                  <template v-slot:activator="{ on }">
                    <v-icon
                      v-on="on"
                      right
                      @click.stop
                      @click="openLink(service.manualLink)"
                      style="
                        display: block;
                        margin-left: auto;
                        margin-right: 10px;
                      "
                      >mdi-help-circle-outline</v-icon
                    >
                  </template>
                  <span>Check out the documentation for more info</span>
                </v-tooltip>
              </div>
            </v-expansion-panel-header>
          </template>
          <span>{{ service.tooltip }}</span>
        </v-tooltip>
        <div
          v-if="
            $store.state.runInfo.active == true &&
            $store.state.runInfo.type == $route.params.workflowName &&
            index == $store.state.runInfo.step
          "
        >
          <v-progress-linear indeterminate color="#1DE9B6"></v-progress-linear>
        </div>
        <v-expansion-panel-content>
          <v-row justify="center">
            <v-col
              v-for="(input, i) in service.Inputs"
              :key="input.name"
              cols="12"
              :xl="input.type === 'combobox' ? 4 : 2"
              :lg="input.type === 'combobox' ? 6 : 3"
              :md="input.type === 'combobox' ? 8 : 4"
              :sm="input.type === 'combobox' ? 12 : 6"
              style="height: fit-content; width: fit-content"
            >
              <v-container :disabled="true" v-if="input.type === 'numeric'"
                ><InputNumeric
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'Inputs'"
              /></v-container>
              <v-container v-if="input.type === 'link'"
                ><InputLink
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'Inputs'"
              /></v-container>
              <v-container v-if="input.type === 'bool'"
                ><InputBool
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'Inputs'"
              /></v-container>
              <v-container v-if="input.type === 'select'"
                ><InputSelect
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'Inputs'"
              /></v-container>
              <v-container v-if="input.type === 'file'"
                ><InputFile
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'Inputs'"
              /></v-container>
              <v-container v-if="input.type === 'boolfile'"
                ><InputBoolFile
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'Inputs'"
              /></v-container>
              <v-container v-if="input.type === 'boolselect'"
                ><InputBoolSelect
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'Inputs'"
              /></v-container>
              <v-container v-if="input.type === 'chip'"
                ><InputChip
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'Inputs'"
              /></v-container>
              <v-container v-if="input.type === 'slide'"
                ><InputSlide
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'Inputs'"
              /></v-container>
              <v-container v-if="input.type === 'combobox'"
                ><InputCombo
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'Inputs'"
              /></v-container>
            </v-col>
          </v-row>
          <v-row v-if="service.extraInputs.length > 0" justify="center">
            <v-btn
              light
              class="mt-5 mb-8"
              style="justify-content: center"
              @click="toggleExtra(index)"
            >
              toggle advance options
            </v-btn></v-row
          >
          <v-divider
            v-if="service.extraInputs.length > 0 && service.showExtra == true"
          ></v-divider>
          <v-row
            style="margin-top: 10px"
            justify="center"
            v-if="service.extraInputs.length > 0 && service.showExtra == true"
          >
            <v-col
              v-for="(input, i) in service.extraInputs"
              :key="input.name"
              cols="12"
              :xl="input.type === 'combobox' ? 4 : 2"
              :lg="input.type === 'combobox' ? 6 : 3"
              :md="input.type === 'combobox' ? 8 : 4"
              :sm="input.type === 'combobox' ? 12 : 6"
              style="height: fit-content; width: fit-content"
            >
              <v-container v-if="input.type === 'numeric'"
                ><InputNumeric
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'extraInputs'"
              /></v-container>
              <v-container v-if="input.type === 'bool'"
                ><InputBool
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'extraInputs'"
              /></v-container>
              <v-container v-if="input.type === 'select'"
                ><InputSelect
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'extraInputs'"
              /></v-container>
              <v-container v-if="input.type === 'file'"
                ><InputFile
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'extraInputs'"
              /></v-container>
              <v-container v-if="input.type === 'boolfile'"
                ><InputBoolFile
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'extraInputs'"
              /></v-container>
              <v-container v-if="input.type === 'boolselect'"
                ><InputBoolSelect
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'extraInputs'"
              /></v-container>
              <v-container v-if="input.type === 'chip'"
                ><InputChip
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'extraInputs'"
              /></v-container>
              <v-container v-if="input.type === 'slide'"
                ><InputSlide
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'extraInputs'"
              /></v-container>
              <v-container v-if="input.type === 'combobox'"
                ><InputCombo
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'extraInputs'"
              /></v-container>
            </v-col>
          </v-row>
        </v-expansion-panel-content>
      </v-expansion-panel>
    </v-expansion-panels>
  </v-row>
</template>

<script>
const { shell } = require("electron");
import InputNumeric from "../components/InputNumeric.vue";
import InputLink from "../components/InputLink.vue";
import InputBool from "../components/InputBool.vue";
import InputSelect from "../components/InputSelect.vue";
import InputFile from "../components/InputFile.vue";
import InputBoolFile from "../components/InputBoolFile.vue";
import InputBoolSelect from "../components/InputBoolSelect.vue";
import InputChip from "../components/InputChip.vue";
import InputSlide from "../components/InputSlide.vue";
import InputCombo from "../components/InputCombo.vue";

export default {
  components: {
    InputLink,
    InputChip,
    InputNumeric,
    InputBool,
    InputSelect,
    InputFile,
    InputBoolFile,
    InputBoolSelect,
    InputSlide,
    InputCombo,
  },
  data: () => ({
    panel: [],
    panel_single: [0],
  }),
  computed: {
    hide() {
      return "display_none";
    },
    services() {
      return this.$store.state[this.$route.params.workflowName];
    },
    inputData() {
      return this.$store.state.data;
    },
  },
  methods: {
    openLink(value) {
      shell.openExternal(value);
    },
    toggleExtra(index) {
      this.$store.commit("toggleExtraCustomWorkflow", {
        workflowName: this.$route.params.workflowName,
        serviceIndex: index,
      });
    },
    check_one(value, index) {
      this.$store.commit("checkCustomService", {
        serviceIndex: index,
        selected: value,
        name: this.$route.params.workflowName,
      });
    },
  },
};
</script>

<style scoped>
.centered-input >>> input {
  text-align: center;
}
.v-card {
  height: 145px;
}
span {
  width: 34px;
}
.container {
  padding: 0;
}
.always {
  padding-left: 58px;
}
.display_none {
  display: none;
}
</style>
