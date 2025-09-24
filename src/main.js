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

Object.defineProperty(Vue.prototype, '$docker', {
  get() {
    return new Docker();
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
    // Gather system specs first, then start Docker monitoring
    this.$store.dispatch('gatherSystemSpecs')
      .then(specs => {
        console.log('System specs gathered:', specs);
        // Start Docker status monitoring after system specs are gathered
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
