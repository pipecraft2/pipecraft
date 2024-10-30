import Vue from "vue"; // Import Vue
import VueRouter from "vue-router";
import { shallowMount, createLocalVue } from "@vue/test-utils";
import Vuex from "vuex";
import Run from "@/components/Run.vue";
import store from "@/store"; // Import the Vuex store
import Vuetify from "vuetify";
// import router from "@/router"; // Import the Vuex store

Vue.use(Vuetify); // Use Vuetify
const localVue = createLocalVue();
localVue.use(Vuex);
localVue.use(VueRouter);

const router = new VueRouter({
  routes: [{ path: "/pipecraft/1", component: Run }],
});

describe("Run.vue", () => {
  it("Runs workflows", () => {
    const wrapper = shallowMount(Run, { store, router, localVue });
    const result = wrapper.vm.toMinsAndSecs(400000);
    expect(result).toBe("6:40");
  });
});
