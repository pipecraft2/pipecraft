<template>
  <draggable
    class="list-group"
    v-model="selectedSteps"
    @start="drag = true"
    @end="push2routeOnEnd($event), (drag = false)"
  >
    <transition-group type="transition" :name="'flip-list'">
      <li
        v-for="(element, index) in selectedSteps"
        class="list-group-item"
        :key="`${element.stepName}${index}`"
      >
        <v-btn
          block
          :style="
            `/step/${element.stepName}/${index}` == $route.path
              ? { borderBottom: 'thin white solid' }
              : { background: '#424242' }
          "
          @click="push2route(element.stepName.replace('/', ''), index)"
        >
          {{ element.stepName }}
          <!-- <v-icon
            v-if="
              $store.state.runInfo.active == false &&
                $store.state.runInfo.type != 'workflow'
            "
            @click.stop="removeAt(index)"
            >mdi-close-box</v-icon
          > -->
          <v-progress-circular
            v-if="
              $store.state.runInfo.active == true &&
              $store.state.runInfo.type == 'workflow' &&
              index == $store.state.runInfo.step
            "
            indeterminate
            color="#1DE9B6"
            :size="20"
          ></v-progress-circular>
          <v-icon v-else @click.stop="removeAt(index)">mdi-close-box</v-icon>
        </v-btn>
      </li>
    </transition-group>
  </draggable>
</template>

<script>
// import About from "../views/About";
import draggable from "vuedraggable";
export default {
  data() {
    return {
      active: "red",
    };
  },

  components: {
    draggable,
  },
  computed: {
    path() {
      return this.$route.path;
    },
    selectedSteps: {
      get() {
        return this.$store.state.selectedSteps;
      },
      set(value) {
        this.$store.commit("DraggableUpdate", value);
      },
    },
  },
  methods: {
    removeAt(index) {
      this.$store.commit("removeStep", index);
      if (this.$route.path != "/home") {
        this.$router.push("/home");
      }
    },
    push2route(stepName, index) {
      let route = `/step/${stepName}/${index}`;
      if (this.$route.path != route) {
        this.$router.push(route);
      }
    },
    push2routeOnEnd(value) {
      let route = `/step/${this.selectedSteps[value.newIndex].stepName}/${
        value.newIndex
      }`;
      if (this.$route.path != route) {
        this.$router.push(route);
      }
    },
  },
};
</script>

<style>
.active {
  color: #1de9b6;
}
,
.flip-list-move {
  transition: transform 0.5s;
}
.no-move {
  transition: transform 0s;
}
.ghost {
  opacity: 0.5;
  background: whitesmoke;
}
.list-group {
  min-height: 20px;
  padding: 0;
}
.list-group-item {
  cursor: pointer;
  padding-top: 10px;
}
.list-group-item i {
  cursor: pointer;
}
li {
  list-style: none;
}
ul.list-group {
  padding-left: 0;
}
.v-btn {
  justify-content: space-between;
}
.v-btn:not(.v-btn--round).v-size--default {
  padding: 0 10px;
}
.routeActive {
  border-bottom: 1px #1de9b6 solid;
  border-bottom-width: 95%;
}
.route {
  background-color: #424242;
}
</style>
