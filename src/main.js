import Vue from "vue";
import App from "./App.vue";
import router from "./router";
import store from "./store";
import './assets/swal.scss';
import vuetify from "./plugins/vuetify";
import { sync } from "vuex-router-sync";
import os from 'os'
sync(store, router);
Vue.config.productionTip = false;

new Vue({
  router,
  store,
  vuetify,
  render: (h) => h(App),
  created() {
    this.$store.commit('setOsType', os.type());
    // Prevent blank screen in Electron builds
    if (this.$route.path != "/home") {
      this.$router.push("/home");
    }
  },
}).$mount("#app");
