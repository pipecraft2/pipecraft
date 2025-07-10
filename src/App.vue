<template>
  <v-app>
    <v-navigation-drawer app clipped width="260" permanent dark>
      <leftNav />
      <!-- <SelectedRoutes /> -->
    </v-navigation-drawer>

    <v-navigation-drawer app clipped right permanent width="72" dark>
      <rightNav />
    </v-navigation-drawer>

    <v-main>
      <transition name="fade">
        <router-view />
      </transition>
    </v-main>

    <DockerPullSnackbar />
  </v-app>
</template>

<script>
import leftNav from "./components/leftNav";
import rightNav from "./components/rightNav";
import DockerPullSnackbar from "./components/DockerPullSnackbar";

export default {
  theme: { dark: true },
  name: "App",

  components: {
    leftNav,
    rightNav,
    DockerPullSnackbar
  },

  data: () => ({
    //
  }),
};
</script>

<style>

::-webkit-scrollbar {
  display: none;
}
.fade-enter-active {
  transition: opacity 1s;
}
.fade-leave-active {
  transition: opacity 0.2;
}
.fade-enter,
.fade-leave-to {
  opacity: 0;
}
#app {
  background: url("assets/bg.jpg") no-repeat center center fixed !important;
  background-size: cover;
}
</style>
