"use strict";

const isDevelopment = process.env.NODE_ENV !== "production";
const log = require("electron-log");
const path = require("path");
import { app, protocol, BrowserWindow, dialog, ipcMain, session } from "electron";
import { createProtocol } from "vue-cli-plugin-electron-builder/lib";
import * as remoteMain from "@electron/remote/main";
import { autoUpdater } from "electron-updater";
remoteMain.initialize();
autoUpdater.logger = require("electron-log");
autoUpdater.logger.transports.file.level = "info";
autoUpdater.autoDownload = false;
var win;
const createDesktopShortcut = require("create-desktop-shortcuts");
const fs = require("fs");
const sudo = require('sudo-prompt');


function sendToRenderer(channel, payload) {
  if (win && !win.isDestroyed()) {
    win.webContents.send(channel, payload);
  }
}

function checkForUpdates() {
  if (isDevelopment) {
    sendToRenderer("update-not-available");
    return Promise.resolve();
  }
  return autoUpdater.checkForUpdates();
}

async function checkScriptsPermissions() {
  const scriptsPath = `${process.resourcesPath}/src/pipecraft-core/service_scripts`;

  if (process.platform === 'linux' && !isDevelopment) {
    try {
      const stats = fs.statSync(scriptsPath);
      const currentPermissions = (stats.mode & 0o777).toString(8);
      
      if (currentPermissions !== '777') {
        await new Promise((resolve, reject) => {
          sudo.exec(
            `chmod -R 777 "${scriptsPath}"`,
            { name: 'PipeCraft' },
            (error, stdout, stderr) => {
              if (error) {
                console.error('Sudo exec error:', error);
                console.error('Stderr:', stderr);
                reject(error);
              } else {
                log.info('Scripts directory permissions updated successfully');
                resolve(stdout);
              }
            }
          );
        });
      }
    } catch (error) {
      log.error('Error checking/setting permissions:', error);
      throw error;
    }
  }
}

log.info("App starting...");

autoUpdater.on("checking-for-update", () => {
  log.info("Checking for update...");
  sendToRenderer("update-checking");
});
autoUpdater.on("update-available", (info) => {
  log.info("Update available:", info && info.version);
  sendToRenderer("update-available", { version: info && info.version });
});
autoUpdater.on("update-not-available", () => {
  log.info("Update not available.");
  sendToRenderer("update-not-available");
});
autoUpdater.on("error", (err) => {
  const message = err && err.message ? err.message : String(err);
  log.error("Error in auto-updater:", message);
  sendToRenderer("update-error", message);
});
autoUpdater.on("download-progress", (progressObj) => {
  log.info("Downloading update:", progressObj.percent);
  sendToRenderer("download-progress", { percent: progressObj.percent });
});
autoUpdater.on("update-downloaded", (info) => {
  log.info("Update downloaded:", info && info.version);
  sendToRenderer("update-downloaded", { version: info && info.version });
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
    minHeight: 700,
    icon: `${process.resourcesPath}/src/pipecraft-core/icon32x32.png`,
    webPreferences: {
      // Use pluginOptions.nodeIntegration, leave this alone
      // See nklayman.github.io/vue-cli-plugin-electron-builder/guide/security.html#node-integration for more info
      nodeIntegration: process.env.ELECTRON_NODE_INTEGRATION,
      contextIsolation: !process.env.ELECTRON_NODE_INTEGRATION,
    },
  });
  remoteMain.enable(win.webContents);
  win.removeMenu();
    // Add keyboard shortcut for DevTools
  win.webContents.on('before-input-event', (event, input) => {
    // Ctrl+Shift+I (Windows/Linux) or Cmd+Opt+I (Mac)
    if ((input.control || input.meta) && input.shift && input.key.toLowerCase() === 'i') {
      win.webContents.toggleDevTools();
    }
  })
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
  try {
    await checkScriptsPermissions();
  } catch (error) {
    log.error('Failed to set permissions:', error);
    // You might want to show a dialog here to inform the user
    dialog.showErrorBox(
      'Permission Error',
      'Failed to set required permissions. Some features might not work correctly.'
    );
  }
  createWindow();
  session.defaultSession.loadExtension(path.join(__dirname, '..', 'devtools5'));
  if (process.env.APPIMAGE && !process.env.WEBPACK_DEV_SERVER_URL) {
    fs.copyFile(
      `${process.resourcesPath}/src/pipecraft-core/icon32x32.png`,
      "/var/tmp/icon32x32.png",
      (err) => {
        if (err) {
          console.log("Error Found:", err);
        }
      }
    );
    createDesktopShortcut({
      onlyCurrentOS: true,
      verbose: true,
      linux: {
        filePath: `${process.env.APPIMAGE}`,
        name: "Pipecraft",
        description: "metabarcoding",
        icon: `/var/tmp/icon32x32.png`,
        type: "Application",
        terminal: false,
        chmod: true,
      },
    });
  }
});

ipcMain.on("update-check", () => {
  checkForUpdates().catch((err) => {
    const message = err && err.message ? err.message : String(err);
    sendToRenderer("update-error", message);
  });
});

ipcMain.on("update-download", () => {
  autoUpdater.downloadUpdate().catch((err) => {
    const message = err && err.message ? err.message : String(err);
    sendToRenderer("update-error", message);
  });
});

ipcMain.on("update-install", () => {
  setImmediate(() => autoUpdater.quitAndInstall());
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
