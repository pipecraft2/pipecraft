<template>
  <div class="text-center">
    <v-menu offset-x :close-on-content-click="false">
      <template v-slot:activator="{ on, attrs }">
        <v-btn
          block
          outlined
          color="white"
          style="background-color: #212121"
          dark
          v-bind="attrs"
          v-on="on"
        >
          add step
        </v-btn>
      </template>
      <v-list dark>
        <v-list-item
          link
          :disabled="Object.values(inputData).includes(item.disabled)"
          v-for="item in items"
          v-bind:key="item.stepName"
          @click="addStep(item, nrOfSelectedSteps)"
        >
          <v-list-item-title>{{ item.stepName }}</v-list-item-title>
          <v-icon>mdi-plus-box</v-icon>
        </v-list-item>
      </v-list>
    </v-menu>
  </div>
</template>

<script>
import { mapState } from "vuex";

export default {
  name: "AddMenu",
  computed: mapState({
    // arrow functions can make the code very succinct!
    inputData: (state) => state.data,
    items: (state) => state.steps,
    nrOfSelectedSteps: (state) => state.selectedSteps.length + 1,
  }),
  methods: {
    addStep(item) {
      this.$store.commit("addStep", {
        step: item,
      });
    },
  },
};
</script>

<style scoped>
.v-btn {
  justify-content: center;
}
.v-list-item__title {
  padding-right: 15px;
  text-transform: uppercase;
}
</style>
