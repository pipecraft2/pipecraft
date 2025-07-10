<template>
  <div class="docker-pull-container">
    <v-snackbar
      v-model="showSnackbar"
      :color="snackbarColor"
      top
      right
      :timeout="-1"
      class="docker-pull-snackbar"
    >
      <div class="d-flex align-center">
        <v-progress-circular
          :value="pullProgress"
          :size="45"
          :width="2"
          color="white"
          class="mr-3"
        >
          {{ pullProgress }}%
        </v-progress-circular>
        <span class="text-body-1">{{ pullStatus }}</span>
      </div>
    </v-snackbar>
  </div>
</template>

<script>
import { mapState } from 'vuex';

export default {
  name: 'DockerPullSnackbar',
  computed: {
    ...mapState({
      pullProgress: state => state.pullProgress,
      pullStatus: state => state.pullStatus,
      pullLoaderActive: state => state.pullLoader.active
    }),
    showSnackbar: {
      get() {
        return this.pullLoaderActive;
      },
      set(value) {
        // Only allow closing if pull is complete or if value is false
        if (this.pullStatus === 'Complete!' || !value) {
          this.$store.commit('deactivatePullLoader');
        }
      }
    },
    snackbarColor() {
      return this.pullStatus === 'Complete!' ? 'primary' : '#363636';
    }
  }
};
</script>

<style scoped>
.docker-pull-container {
  position: fixed;
  top: 20px;
  right: 72px;
  z-index: 9999;
}

.test-button {
  margin: 0px;
}

.docker-pull-snackbar {
  z-index: 9999 !important;
}

::v-deep .v-snack__wrapper {
  box-shadow: none !important;
  border-bottom: 1px solid grey !important;
  border-radius: 0px !important;
  right: 64px !important;
  height: 84px !important;
  margin-top: 0 !important;
  min-width: 250px !important;
}
</style> 