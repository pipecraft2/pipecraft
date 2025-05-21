<template>
  <v-container fluid class="text-center">
    <v-row justify="center" align="center">
      <v-col cols="12" class="d-flex flex-column align-center">
        <!-- Circular Progress with Percentage -->
        <v-progress-circular
          :size="70"
          :width="7"
          :value="pullProgress"
          color="purple"
          class="mb-4"
        >
          {{ pullProgress }}%
        </v-progress-circular>

        <!-- Status Text -->
        <div class="text-subtitle-1">
          {{ pullStatus || 'Preparing download...' }}
        </div>
      </v-col>
    </v-row>
  </v-container>
</template>

<script>
import { mapState } from 'vuex';

export default {
  name: 'PullLoader',
  computed: {
    ...mapState({
      pullProgress: state => state.pullProgress,
      pullStatus: state => state.pullStatus
    })
  }
};
</script>

<style scoped>
.v-progress-circular {
  margin-bottom: 16px;
}
</style>
