<template>
  <v-app>
    <v-navigation-drawer app clipped width="260" permanent dark>
      <leftNav />
      <!-- <SelectedRoutes /> -->
    </v-navigation-drawer>

    <v-navigation-drawer app clipped right permanent width="72" dark>
      <rightNav />
      <template #append>
        <UpdateButton />
      </template>
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
import UpdateButton from "./components/UpdateButton";
import DockerPullSnackbar from "./components/DockerPullSnackbar";

export default {
  theme: { dark: true },
  name: "App",

  components: {
    leftNav,
    rightNav,
    UpdateButton,
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
html,
body,
#app,
.v-application,
.v-application--wrap {
  min-height: 100vh;
  height: 100%;
}

.v-application,
.v-application--wrap {
  background: url("assets/bg.jpg") no-repeat center center !important;
  background-size: cover !important;
  background-color: #000 !important;
}
</style>
