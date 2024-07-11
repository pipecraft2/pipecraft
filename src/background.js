"use strict";

const isDevelopment = process.env.NODE_ENV !== "production";
const log = require("electron-log");
import { app, protocol, BrowserWindow, dialog, ipcMain } from "electron";
import { createProtocol } from "vue-cli-plugin-electron-builder/lib";
import installExtension, { VUEJS_DEVTOOLS } from "electron-devtools-installer";
import * as remoteMain from "@electron/remote/main";
import { autoUpdater } from "electron-updater";
remoteMain.initialize();
autoUpdater.logger = require("electron-log");
autoUpdater.logger.transports.file.level = "info";
autoUpdater.autoDownload = false;
var win;
const createDesktopShortcut = require('create-desktop-shortcuts');
const fs = require("fs");


function update() {
  autoUpdater.checkForUpdates();
}

log.info("App starting...");

autoUpdater.on("checking-for-update", () => {
  return log.info("Checking for update...");
});
autoUpdater.on("update-available", () => {
  dialog
    .showMessageBox({
      type: "question",
      title: "Found Updates",
      message: "A newer version of pipecraft is available",
      buttons: ["Update now", "Cancel"],
    })
    .then((button) => {
      console.log(button);
      if (button.response === 0) {
        console.log(button);
        autoUpdater.downloadUpdate();
      }
    });
  return log.info("Update available.");
});
autoUpdater.on("update-not-available", () => {
  dialog.showMessageBox({
    type: "question",
    title: "No updates found",
    message: "Pipecraft is up to date",
    buttons: ["Ok"],
  });
  win.webContents.send("update-not-available");
  return log.info("Update not available.");
});
autoUpdater.on("error", (err) => {
  console.log(err);
  dialog.showMessageBox({
    type: "error",
    title: "Something went wrong",
    message: `${err}`,
  });
  win.webContents.send("update-error", err);
  return log.info("Error in auto-updater. " + err);
});
autoUpdater.on("download-progress", (progressObj) => {
  console.log(progressObj);
  return log.info("downloading update");
});
autoUpdater.on("update-downloaded", () => {
  win.webContents.send("update-downloaded");
  log.info("Update downloaded");
  dialog
    .showMessageBox({
      title: "Install Updates",
      message: "Updates downloaded, Pipecraft will restart to install updates.",
    })
    .then(() => {
      setImmediate(() => autoUpdater.quitAndInstall());
    });
});

// Scheme must be registered before the app is ready
protocol.registerSchemesAsPrivileged([
  { scheme: "app", privileges: { secure: true, standard: true } },
]);

async function createWindow() {
  // Create the browser window.
  win = new BrowserWindow({
    width: 1200,
    height: 700,
    minWidth: 1200,
    icon: `${process.resourcesPath}/src/pipecraft-core/icon32x32.png`,
    webPreferences: {
      // Use pluginOptions.nodeIntegration, leave this alone
      // See nklayman.github.io/vue-cli-plugin-electron-builder/guide/security.html#node-integration for more info
      nodeIntegration: process.env.ELECTRON_NODE_INTEGRATION,
      contextIsolation: !process.env.ELECTRON_NODE_INTEGRATION,
      enableRemoteModule: true,
    },
  });
  win.webContents.openDevTools();
  win.removeMenu();
  win.webContents.on("new-window", function (e, url) {
    e.preventDefault();
    require("electron").shell.openExternal(url);
  });

  if (process.env.WEBPACK_DEV_SERVER_URL) {
    // Load the url of the dev server if in development mode
    await win.loadURL(process.env.WEBPACK_DEV_SERVER_URL);
    if (!process.env.IS_TEST) win.webContents.openDevTools();
  } else {
    createProtocol("app");
    // Load the index.html when not in development
    win.loadURL("app://./index.html");
  }
}

// Quit when all windows are closed.
app.on("window-all-closed", () => {
  // On macOS it is common for applications and their menu bar
  // to stay active until the user quits explicitly with Cmd + Q
  if (process.platform !== "darwin") {
    app.quit();
  }
});

app.on("activate", () => {
  // On macOS it's common to re-create a window in the app when the
  // dock icon is clicked and there are no other windows open.
  if (BrowserWindow.getAllWindows().length === 0) createWindow();
});

// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on("ready", async () => {
  if (isDevelopment && !process.env.IS_TEST) {
    // Install Vue Devtools
    try {
      await installExtension("iaajmlceplecbljialhhkmedjlpdblhp");
    } catch (e) {
      console.error("Vue Devtools failed to install:", e.toString());
    }
  }
  createWindow();
  if (isDevelopment != true) {
    if (process.env.APPIMAGE.includes('AppImage')) {
      fs.copyFile(`${process.resourcesPath}/src/pipecraft-core/icon32x32.png`, "/var/tmp/icon32x32.png", (err) => {
        if (err) {
          console.log("Error Found:", err);
        }
      });
      createDesktopShortcut({
        onlyCurrentOS: true,
        verbose: true,
        linux: {
          filePath: `${process.env.APPIMAGE}`,
          name: 'Pipecraft',
          description: 'My app description',
          icon: `/var/tmp/icon32x32.png`,
          type: 'Application',
          terminal: false,
          chmod: true,
        },
      })
    }
  }
});

ipcMain.on("update", () => {
  update();
  console.log("updating");
});

// Exit cleanly on request from parent process in development mode.
if (isDevelopment) {
  if (process.platform === "win32") {
    process.on("message", (data) => {
      if (data === "graceful-exit") {
        app.quit();
      }
    });
  } else {
    process.on("SIGTERM", () => {
      app.quit();
    });
  }
}
