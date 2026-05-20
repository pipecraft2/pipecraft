<template>
  <v-menu
    v-model="popupOpen"
    :close-on-content-click="false"
    left
    offset-x
    :nudge-left="10"
    transition="slide-x-reverse-transition"
  >
    <template v-slot:activator="{ on: menuOn, attrs }">
      <v-list-item ripple link v-bind="attrs" style="margin-inline: 8px;">
        <v-tooltip left nudge-left="10" :disabled="popupOpen">
          <template v-slot:activator="{ on: tipOn }">
            <v-list-item-content
              v-on="tipOn"
              @click="onClick($event, menuOn)"
            >
              <v-badge
                :value="showBadge"
                color="#1DE9B6"
                dot
                overlap
                class="update-badge"
              >
                <v-progress-circular
                  v-if="state === 'checking' || state === 'downloading'"
                  indeterminate
                  color="#1DE9B6"
                  size="24"
                />
                <v-icon v-else :color="iconColor">mdi-update</v-icon>
              </v-badge>
            </v-list-item-content>
          </template>
          <span>{{ tooltip }}</span>
        </v-tooltip>
      </v-list-item>
    </template>

    <v-card color="#212121" dark width="280" elevation="8">
      <v-card-title class="py-2 px-3">
        <span class="update-popup-title">{{ popupTitle }}</span>
        <v-spacer />
        <v-btn icon x-small @click="closePopup">
          <v-icon small>mdi-close</v-icon>
        </v-btn>
      </v-card-title>

      <v-card-text class="py-2 px-3 update-popup-body">
        <div>Current: <b>v{{ currentVersion }}</b></div>
        <div v-if="newVersion">
          Latest: <b class="update-latest">v{{ newVersion }}</b>
        </div>
        <div v-if="state === 'error'" class="update-error">{{ error }}</div>
        <div v-if="state === 'downloading' && downloadPercent != null" class="mt-2">
          <v-progress-linear
            :value="downloadPercent"
            color="#1DE9B6"
            height="6"
            rounded
          />
          <div class="mt-1">{{ Math.round(downloadPercent) }}% downloaded</div>
        </div>
        <a
          v-if="newVersion && state !== 'up-to-date'"
          href="#"
          class="update-notes-link"
          @click.prevent="openReleaseNotes"
        >Release notes</a>
      </v-card-text>

      <v-card-actions v-if="state === 'available'" class="px-3 pb-3">
        <v-btn x-small text @click="skipVersion">Skip</v-btn>
        <v-spacer />
        <v-btn x-small text @click="closePopup">Later</v-btn>
        <v-btn x-small color="#1DE9B6" dark @click="install">Install</v-btn>
      </v-card-actions>
      <v-card-actions v-else-if="state === 'ready'" class="px-3 pb-3">
        <v-spacer />
        <v-btn x-small color="#1DE9B6" dark @click="restartNow">Restart now</v-btn>
      </v-card-actions>
      <v-card-actions
        v-else-if="state === 'up-to-date' || state === 'error'"
        class="px-3 pb-3"
      >
        <v-spacer />
        <v-btn x-small text @click="closePopup">OK</v-btn>
      </v-card-actions>
    </v-card>
  </v-menu>
</template>

<script>
const { ipcRenderer, shell } = require("electron");
const { app } = require("@electron/remote");

const SKIP_KEY = "pipecraft.update.skippedVersion";
const LAUNCH_CHECK_DELAY_MS = 4000;

export default {
  name: "UpdateButton",
  data() {
    return {
      state: "idle",
      popupOpen: false,
      currentVersion: app.getVersion(),
      newVersion: null,
      error: null,
      downloadPercent: null,
      silentLaunchCheck: true,
    };
  },
  computed: {
    showBadge() {
      return this.state === "available" || this.state === "ready";
    },
    iconColor() {
      if (this.state === "available" || this.state === "ready") return "#1DE9B6";
      if (this.state === "error") return "#FF7043";
      return "white";
    },
    tooltip() {
      switch (this.state) {
        case "checking":
          return "Checking for updates...";
        case "available":
          return `Update available: v${this.newVersion}`;
        case "downloading":
          return "Downloading update...";
        case "ready":
          return "Update ready — restart to install";
        case "up-to-date":
          return "PipeCraft is up to date";
        case "error":
          return "Update check failed";
        default:
          return "Check for updates";
      }
    },
    popupTitle() {
      switch (this.state) {
        case "available":
          return "Update available";
        case "ready":
          return "Update ready";
        case "up-to-date":
          return "Up to date";
        case "error":
          return "Update check failed";
        case "downloading":
          return "Downloading update";
        default:
          return "Updates";
      }
    },
  },
  mounted() {
    ipcRenderer.on("update-checking", this.onChecking);
    ipcRenderer.on("update-available", this.onAvailable);
    ipcRenderer.on("update-not-available", this.onNotAvailable);
    ipcRenderer.on("download-progress", this.onDownloadProgress);
    ipcRenderer.on("update-downloaded", this.onDownloaded);
    ipcRenderer.on("update-error", this.onError);

    setTimeout(() => {
      ipcRenderer.send("update-check");
    }, LAUNCH_CHECK_DELAY_MS);
  },
  beforeDestroy() {
    ipcRenderer.removeListener("update-checking", this.onChecking);
    ipcRenderer.removeListener("update-available", this.onAvailable);
    ipcRenderer.removeListener("update-not-available", this.onNotAvailable);
    ipcRenderer.removeListener("download-progress", this.onDownloadProgress);
    ipcRenderer.removeListener("update-downloaded", this.onDownloaded);
    ipcRenderer.removeListener("update-error", this.onError);
  },
  methods: {
    onChecking() {
      this.state = "checking";
      this.error = null;
    },
    onAvailable(_event, info) {
      this.newVersion = info && info.version ? info.version : null;
      if (this.newVersion && localStorage.getItem(SKIP_KEY) === this.newVersion) {
        this.state = "idle";
        return;
      }
      this.state = "available";
      this.popupOpen = true;
    },
    onNotAvailable() {
      this.state = "up-to-date";
      this.newVersion = null;
      if (this.silentLaunchCheck) return;
      this.popupOpen = true;
    },
    onDownloadProgress(_event, progress) {
      this.state = "downloading";
      this.downloadPercent =
        progress && progress.percent != null ? progress.percent : null;
      this.popupOpen = true;
    },
    onDownloaded(_event, info) {
      this.newVersion = (info && info.version) || this.newVersion;
      this.state = "ready";
      this.downloadPercent = 100;
      this.popupOpen = true;
    },
    onError(_event, message) {
      this.state = "error";
      this.error = message || "Unknown error";
      if (this.silentLaunchCheck) return;
      this.popupOpen = true;
    },
    onClick(event, menuOn) {
      const idle = this.state === "idle";
      if (idle) {
        this.silentLaunchCheck = false;
        this.state = "checking";
        this.error = null;
        ipcRenderer.send("update-check");
        return;
      }
      // Forward the click to v-menu's activator so it toggles
      if (menuOn && typeof menuOn.click === "function") {
        menuOn.click(event);
      } else {
        this.popupOpen = !this.popupOpen;
      }
    },
    install() {
      this.state = "downloading";
      this.downloadPercent = 0;
      ipcRenderer.send("update-download");
    },
    restartNow() {
      ipcRenderer.send("update-install");
    },
    skipVersion() {
      if (this.newVersion) {
        localStorage.setItem(SKIP_KEY, this.newVersion);
      }
      this.popupOpen = false;
      this.state = "idle";
    },
    closePopup() {
      this.popupOpen = false;
      if (this.state === "up-to-date") {
        this.state = "idle";
      }
    },
    openReleaseNotes() {
      if (!this.newVersion) return;
      shell.openExternal(
        `https://github.com/pipecraft2/pipecraft/releases/tag/v${this.newVersion}`
      );
    },
  },
};
</script>

<style scoped>
.update-badge {
  justify-content: center;
  width: 100%;
}
.update-popup-title {
  font-size: 14px;
}
.update-popup-body {
  font-size: 12px;
  line-height: 1.5;
}
.update-latest {
  color: #1de9b6;
}
.update-error {
  color: #ff7043;
  margin-top: 4px;
}
.update-notes-link {
  color: #80cbc4;
  font-size: 11px;
  display: inline-block;
  margin-top: 6px;
  text-decoration: none;
}
.update-notes-link:hover {
  text-decoration: underline;
}
</style>
