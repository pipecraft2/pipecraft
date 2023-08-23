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
    <v-card-actions style="justify-content: center">
      <v-row style="justify-content: center">
        <v-tooltip top>
          <template v-slot:activator="{ on }">
            <v-col v-on="on" style="padding: 0" cols="10" offset="0">
              <v-select
                @change="inputUpdate(input.value)"
                v-model="input.value"
                :items="input.items"
                :label="input.name"
                multiple
              >
                <template v-slot:prepend-item>
                  <v-list-item ripple @click="toggle">
                    <v-list-item-action>
                      <v-icon
                        :color="input.value.length > 0 ? 'indigo darken-4' : ''"
                      >
                        {{ icon }}
                      </v-icon>
                    </v-list-item-action>
                    <v-list-item-content>
                      <v-list-item-title> Select All </v-list-item-title>
                    </v-list-item-content>
                  </v-list-item>
                  <v-divider class="mt-2"></v-divider>
                </template>

                <template v-slot:selection="{ item, index }">
                  <v-chip v-if="index === 0">
                    <span>{{ item }}</span>
                  </v-chip>
                  <v-chip v-if="index === 1">
                    <span>{{ item }}</span>
                  </v-chip>
                  <v-chip v-if="index === 2">
                    <span>{{ item }}</span>
                  </v-chip>
                  <v-chip v-if="index === 3">
                    <span>{{ item }}</span>
                  </v-chip>
                  <v-chip v-if="index === 4">
                    <span>{{ item }}</span>
                  </v-chip>
                  <v-chip v-if="index === 5">
                    <span>{{ item }}</span>
                  </v-chip>
                  <v-chip v-if="index === 6">
                    <span>{{ item }}</span>
                  </v-chip>
                  <span v-if="index === 6" class="grey--text text-caption">
                    (+{{ input.value.length - 7 }} others)
                  </span>
                </template>
              </v-select>
            </v-col>
          </template>
          <span>{{ input.tooltip }}</span>
        </v-tooltip>
      </v-row>
    </v-card-actions>
  </v-card>
</template>

<script>
export default {
  data: () => ({}),

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
    All() {
      return this.input.value.length === this.input.items.length;
    },
    Some() {
      return this.input.value.length > 0 && !this.All;
    },
    icon() {
      if (this.All) return "mdi-close-box";
      if (this.Some) return "mdi-minus-box";
      return "mdi-checkbox-blank-outline";
    },
  },
  methods: {
    toggle() {
      this.$nextTick(() => {
        if (this.All) {
          this.input.value = [];
          this.inputUpdate(this.input.value);
        } else {
          this.input.value = this.input.items.slice();
          this.inputUpdate(this.input.value);
        }
      });
    },
    inputUpdate(value) {
      if (this.$route.params.workflowName) {
        console.log(value);
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

<style scoped>
.v-text-field input {
  text-align: center;
}
div >>> div.v-select__selections {
  margin-top: 8px;
}
.v-text-field {
  padding-top: 20px;
}
</style>
