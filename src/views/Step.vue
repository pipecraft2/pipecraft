<template>
  <!-- each sub component must recive -->
  <v-row justify="center" style="padding-top: 20px">
    <v-expansion-panels dark multiple popout>
      <v-expansion-panel v-for="(service, index) in services" :key="index">
        <v-tooltip top>
          <template v-slot:activator="{ on }">
            <v-expansion-panel-header v-on="on" style="padding-left: 40%">
              <v-checkbox
                @change="check_one($event, index)"
                @click.stop
                v-model="service.selected"
                style="max-width: 34px; padding-right: 10px"
              ></v-checkbox
              >{{ service.serviceName.toUpperCase().replaceAll('_', ' ') }}
            </v-expansion-panel-header>
          </template>
          <span>{{ service.tooltip }}</span>
        </v-tooltip>
        <v-expansion-panel-content>
          <v-row justify="center">
            <v-col
              v-for="(input, i) in service.Inputs"
              :key="input.name"
              cols="12"
              :xl="input.type === 'combobox' ? 4 : 2"
              :lg="input.type === 'combobox' ? 6 : 3"
              :md="input.type === 'combobox' ? 8 : 4"
              :sm="input.type === 'combobox' ? 12 : 3"
              style="height: fit-content; width: fit-content"
            >
              <v-container v-if="input.type === 'combobox'"
                ><InputCombo
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'Inputs'"
              /></v-container>
              <v-container v-if="input.type === 'numeric'"
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
            </v-col>
          </v-row>
          <v-row v-if="service.extraInputs.length > 0" justify="center">
            <v-btn
              light
              class="mt-5 mb-8"
              style="justify-content: center"
              @click="toggleExtra(index)"
            >
              toggle advanced options
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
              :sm="input.type === 'combobox' ? 12 : 3"
              style="height: fit-content; width: fit-content"
            >
              <v-container v-if="input.type === 'combobox'"
                ><InputCombo
                  :serviceIndex="index"
                  :inputIndex="i"
                  :list="'extraInputs'"
              /></v-container>
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
            </v-col>
          </v-row>
        </v-expansion-panel-content>
      </v-expansion-panel>
    </v-expansion-panels>
  </v-row>
</template>

<script>
import InputNumeric from "../components/InputNumeric.vue";
import InputBool from "../components/InputBool.vue";
import InputSelect from "../components/InputSelect.vue";
import InputFile from "../components/InputFile.vue";
import InputBoolFile from "../components/InputBoolFile.vue";
import InputBoolSelect from "../components/InputBoolSelect.vue";
import InputChip from "../components/InputChip.vue";
import InputSlide from "../components/InputSlide.vue";
import InputCombo from "../components/InputCombo.vue";
import InputLink from "../components/InputLink.vue";

export default {
  name: "Home",
  components: {
    InputCombo,
    InputChip,
    InputNumeric,
    InputBool,
    InputSelect,
    InputFile,
    InputBoolFile,
    InputBoolSelect,
    InputSlide,
    InputLink,
  },
  computed: {
    services() {
      return this.$store.state.selectedSteps[this.$route.params.order].services;
    },
  },
  methods: {
    toggleExtra(index) {
      this.$store.commit("toggleExtra", {
        stepIndex: this.$route.params.order,
        serviceIndex: index,
      });
    },
    check_one(value, index) {
      this.$store.commit("checkService", {
        stepIndex: this.$route.params.order,
        serviceIndex: index,
        selected: value,
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
</style>
