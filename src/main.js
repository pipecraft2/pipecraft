import Vue from "vue";
import App from "./App.vue";
import router from "./router";
import store from "./store";
import './assets/swal.scss';
import vuetify from "./plugins/vuetify";
import { sync } from "vuex-router-sync";
import os from 'os'
import '@mdi/font/css/materialdesignicons.css';
const Docker = require('dockerode');
const { getDockerodeOptionsFromContextSync } = require("./utils/dockerContext");

let dockerOptions;

try {
  dockerOptions = getDockerodeOptionsFromContextSync();
} catch (error) {
  console.error(error.message);
  dockerOptions = {};
}

Object.defineProperty(Vue.prototype, '$docker', {
  get() {
    return new Docker(dockerOptions);
  }
});

sync(store, router);
Vue.config.productionTip = false;

new Vue({
  router,
  store,
  vuetify,
  render: (h) => h(App),
  created() {
    // Gather system specs first, then start Docker monitoring.
    // fetchDockerInfo is dispatched by the monitor when Docker becomes
    // "running" (and again when opening Resource Manager), so we don't
    // race Docker Desktop boot with an eager info fetch here.
    this.$store.dispatch('gatherSystemSpecs')
      .then(specs => {
        console.log('System specs gathered:', specs);
        this.$store.dispatch('startDockerStatusMonitoring');
      })
      .catch(error => {
        console.error('Failed to gather system specs:', error);
        // Still try to start Docker monitoring even if system specs fail
        this.$store.dispatch('startDockerStatusMonitoring');
      });
    this.$store.commit('setOsType', os.type());
    // Prevent blank screen in Electron builds
    if (this.$route.path != "/home") {
      this.$router.push("/home");
    }
  },
}).$mount("#app");
