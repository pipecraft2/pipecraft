<template>
  <v-dialog
    :value="visible"
    persistent
    max-width="560"
    overlay-opacity="0.7"
  >
    <v-card dark>
      <v-card-title class="text-h6">
        {{ dialogTitle }}
      </v-card-title>
      <v-card-text>
        <div v-if="starting" class="d-flex align-center">
          <v-progress-circular
            indeterminate
            color="primary"
            size="28"
            width="3"
            class="mr-4"
            aria-label="Starting container engine"
          ></v-progress-circular>
          <span>{{ startingText }}</span>
        </div>
        <template v-else-if="needInstallWarning">
          <v-alert type="warning" dense class="mb-4">
            PipeCraft runs every analysis step in a container. Without Docker or
            Podman you can still browse the interface, but you cannot start a
            workflow.
          </v-alert>
          <p class="mb-3">
            Install one engine, then restart PipeCraft or click Check again.
          </p>
          <div class="d-flex flex-column">
            <v-btn
              text
              color="primary"
              class="justify-start px-0"
              @click="handleOpenLink(dockerInstallUrl)"
            >
              Docker installation
              <v-icon right small>mdi-open-in-new</v-icon>
            </v-btn>
            <v-btn
              text
              color="primary"
              class="justify-start px-0"
              @click="handleOpenLink(podmanInstallUrl)"
            >
              Podman installation
              <v-icon right small>mdi-open-in-new</v-icon>
            </v-btn>
          </div>
        </template>
        <template v-else-if="needChooser">
          <p class="mb-4">
            PipeCraft found both Docker and Podman. Only one engine should be
            used. This choice can be changed later in Resource Manager.
          </p>
          <v-alert
            v-if="startError"
            type="error"
            dense
            class="mb-4"
          >
            {{ startError }}
          </v-alert>
          <v-row>
            <v-col cols="12" sm="6">
              <div class="engine-option">
                <div class="text-subtitle-1">Docker</div>
                <div class="engine-status">{{ dockerStatusLabel }}</div>
                <v-btn
                  class="mt-3"
                  color="primary"
                  :disabled="starting || !dockerInstalled"
                  @click="handleChoose('docker')"
                >
                  Use Docker
                </v-btn>
              </div>
            </v-col>
            <v-col cols="12" sm="6">
              <div class="engine-option">
                <div class="text-subtitle-1">Podman</div>
                <div class="engine-status">{{ podmanStatusLabel }}</div>
                <v-btn
                  class="mt-3"
                  color="primary"
                  :disabled="starting || !podmanInstalled"
                  @click="handleChoose('podman')"
                >
                  Use Podman
                </v-btn>
              </div>
            </v-col>
          </v-row>
          <v-checkbox
            v-model="rememberAsDefault"
            dark
            hide-details
            class="mt-4"
            label="Use as default"
          ></v-checkbox>
          <div class="caption mt-1" style="opacity: 0.75">
            When checked, this engine starts automatically next time. Uncheck to
            choose again at launch.
          </div>
        </template>
        <v-alert v-else-if="startError" type="error" dense>
          {{ startError }}
        </v-alert>
      </v-card-text>
      <v-card-actions v-if="needInstallWarning && !starting">
        <v-spacer></v-spacer>
        <v-btn text @click="handleDismissInstallWarning">Continue anyway</v-btn>
        <v-btn color="primary" text @click="handleCheckAgain">Check again</v-btn>
      </v-card-actions>
      <v-card-actions v-else-if="startError && !starting && !needChooser">
        <v-spacer></v-spacer>
        <v-btn text color="primary" @click="handleRetry">Retry</v-btn>
      </v-card-actions>
    </v-card>
  </v-dialog>
</template>

<script>
const { shell } = require("electron");
import { mapState } from "vuex";

export default {
  name: "ContainerEngineDialog",
  data() {
    return {
      rememberAsDefault: true,
      dockerInstallUrl: "https://docs.docker.com/get-docker/",
      podmanInstallUrl: "https://podman.io/docs/installation",
    };
  },
  computed: {
    ...mapState({
      needChooser: (state) => state.containerRuntime.needChooser,
      needInstallWarning: (state) => state.containerRuntime.needInstallWarning,
      starting: (state) => state.containerRuntime.starting,
      startError: (state) => state.containerRuntime.startError,
      sessionEngine: (state) => state.containerRuntime.sessionEngine,
      dockerInstalled: (state) => state.containerRuntime.available.docker,
      podmanInstalled: (state) => state.containerRuntime.available.podman,
      dockerRunning: (state) => state.containerRuntime.running.docker,
      podmanRunning: (state) => state.containerRuntime.running.podman,
    }),
    visible() {
      return (
        this.needChooser ||
        this.needInstallWarning ||
        this.starting ||
        Boolean(this.startError)
      );
    },
    dialogTitle() {
      if (this.starting) {
        return "Starting container engine";
      }
      if (this.needInstallWarning) {
        return "No container engine found";
      }
      if (this.needChooser) {
        return "Choose a container engine";
      }
      return "Container engine";
    },
    startingText() {
      const name = this.sessionEngine === "podman" ? "Podman" : "Docker";
      return `Starting ${name}. This can take up to a minute.`;
    },
    dockerStatusLabel() {
      if (!this.dockerInstalled) {
        return "Not installed";
      }
      return this.dockerRunning ? "Already running" : "Installed, currently stopped";
    },
    podmanStatusLabel() {
      if (!this.podmanInstalled) {
        return "Not installed";
      }
      return this.podmanRunning ? "Already running" : "Installed, currently stopped";
    },
  },
  methods: {
    handleOpenLink(url) {
      shell.openExternal(url);
    },
    handleDismissInstallWarning() {
      this.$store.commit("setNeedInstallWarning", false);
    },
    async handleCheckAgain() {
      try {
        await this.$store.dispatch("bootstrapContainerRuntime");
      } catch (error) {
        console.error("Failed to start container engine:", error);
      }
    },
    async handleChoose(engine) {
      try {
        await this.$store.dispatch("activateEngine", {
          engine,
          persistDefault: this.rememberAsDefault,
          start: true,
        });
      } catch (error) {
        console.error("Failed to start container engine:", error);
      }
    },
    async handleRetry() {
      if (!this.sessionEngine) {
        return;
      }
      try {
        await this.$store.dispatch("activateEngine", {
          engine: this.sessionEngine,
          persistDefault: this.rememberAsDefault,
          start: true,
        });
      } catch (error) {
        console.error("Failed to start container engine:", error);
      }
    },
  },
};
</script>

<style scoped>
.engine-option {
  border: 1px solid rgba(255, 255, 255, 0.2);
  padding: 16px;
  min-height: 150px;
}
.engine-status {
  opacity: 0.75;
  font-size: 13px;
  margin-top: 4px;
}
</style>
