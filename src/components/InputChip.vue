<template>
  <v-tooltip top>
    <template v-slot:activator="{ on }">
      <div v-on="on">
        <v-card
          :disabled="
            Object.values(inputData).includes(input.disabled) ||
            $store.state.runInfo.active == true
          "
          light
          elevation="2"
          style="height: fit-content; resize: auto; min-height: 145px"
        >
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

          <v-card-actions style="justify-content: center">
            <v-row style="height: fit-content; resize: auto"
              ><v-col style="padding: 0" cols="10" offset="1">
                <v-combobox
                  type="text"
                  @change="inputUpdate(input.value)"
                  append-icon=""
                  deletable-chips
                  small-chips
                  multiple
                  :rules="input.rules"
                  style="text-transform: uppercase"
                  @keydown="input.iupac ? IUPAC($event) : void 0"
                  v-model="input.value"
                ></v-combobox>
              </v-col>
            </v-row>
          </v-card-actions>
        </v-card>
      </div>
    </template>
    <span>{{ input.tooltip }}</span>
  </v-tooltip>
</template>

<script>
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
    IUPAC(event) {
      if (
        ![
          "backspace",
          "i",
          "a",
          "c",
          "g",
          "t",
          "r",
          "y",
          "s",
          "w",
          "k",
          "m",
          "b",
          "d",
          "h",
          "v",
          "n",
        ].includes(event.key.toLowerCase())
      ) {
        event.preventDefault();
      }
    },
    inputUpdate(value) {
      if (this.$route.params.workflowName) {
        this.$store.commit("premadeInputUpdate", {
          workflowName: this.$route.params.workflowName,
          serviceIndex: this.$attrs.serviceIndex,
          inputIndex: this.$attrs.inputIndex,
          listName: this.$attrs.list,
          value: value,
        });
      } else {
        this.$store.commit("inputUpdate", {
          stepIndex: this.$route.params.order,
          serviceIndex: this.$attrs.serviceIndex,
          inputIndex: this.$attrs.inputIndex,
          listName: this.$attrs.list,
          value: value,
        });
      }
    },
  },
};
</script>

<style>
.v-text-field input {
  text-align: center;
}
.v-text-field__details {
  position: relative;
  margin-bottom: 5%;
}
</style>
