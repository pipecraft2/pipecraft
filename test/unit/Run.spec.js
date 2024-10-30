import { shallowMount, createLocalVue } from "@vue/test-utils";
import Vuex from "vuex";
import Run from "@/components/Run.vue";
import store from "@/store"; // Import the Vuex store

const localVue = createLocalVue();
localVue.use(Vuex);

describe("Run.vue", () => {
  it("Runs workflows", () => {
    const wrapper = shallowMount(Run, { store, localVue });
    expect(wrapper.text()).toBe("Hello from Vuex store!");

    const result = wrapper.vm.toMinsAndSecs(400000);
    expect(result).toBe("6:40");
  });
});
