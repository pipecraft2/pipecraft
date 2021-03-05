"use strict";

import { app, protocol, BrowserWindow, ipcMain } from "electron";
import { createProtocol } from "vue-cli-plugin-electron-builder/lib";
import installExtension, { VUEJS_DEVTOOLS } from "electron-devtools-installer";
import { promises as fs } from "fs";
// import { exitCode, stdout } from "process";
var Docker = require("dockerode");
var docker = new Docker({ socketPath: "//./pipe/docker_engine" });
const streams = require("memory-streams");
var stdout = new streams.WritableStream();
var stderr = new streams.WritableStream();

const isDevelopment = process.env.NODE_ENV !== "production";

function runStep(serviceImage, scriptName, envVariables) {
  return new Promise((resolve, reject) => {
    // var imageList = docker.listImages();
    // resolve(image)
    var result = docker
      .run("pipecraft/mothur:1.43", ["bash", "-c", "ls"], [stdout, stderr], {
        Tty: false,
      })
      .then(function(data) {
        console.log("stdout:", stdout.toString());
        console.log("stderr:", stderr.toString());
        var output = data[0];
        var container = data[1];
        console.log(output);
        container.remove();
        return "tere";
      })
      .then(function(data) {
        console.log("container removed");
      })
      .catch(function(err) {
        console.log(err);
      });
    resolve(result);
  });
}

ipcMain.on(
  "runStep",
  async (event, imageName, scriptName, envVariables, Input) => {
    var result = await docker
      .run(
        imageName,
        ["ash", "-c", `ls && /scripts/${scriptName}`],
        [stdout, stderr],
        {
          Tty: false,
          WorkingDir: "/input",
          Volumes: {},
          HostConfig: {
            Binds: [
              "C:\\Users\\m_4_r\\Desktop\\pipecraft-vue\\pipecraft-core\\service_scripts:/scripts", // Edit path for build
              `${Input}:/input`,
            ],
          },
          Env: envVariables,
        },
      )
      .then(([res, container]) => {
        console.log(res);
        console.log("stdout: %j", stdout.toString());
        console.log("stderr: %j", stderr.toString());
        container.remove();
        if (res.StatusCode === 0) {
          return stdout.toString();
        } else {
          return stderr.toString();
        }
      })
      .catch((err) => {
        console.log(err);
        return err;
      });
    event.returnValue = result;
    stdout = new streams.WritableStream();
    stderr = new streams.WritableStream();
  },
);

// Scheme must be registered before the app is ready
protocol.registerSchemesAsPrivileged([
  { scheme: "app", privileges: { secure: true, standard: true } },
]);

async function createWindow() {
  // Create the browser window.
  const win = new BrowserWindow({
    width: 800,
    height: 600,
    webPreferences: {
      // Use pluginOptions.nodeIntegration, leave this alone
      // See nklayman.github.io/vue-cli-plugin-electron-builder/guide/security.html#node-integration for more info
      nodeIntegration: process.env.ELECTRON_NODE_INTEGRATION,
    },
  });
  // win.removeMenu();

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
      await installExtension(VUEJS_DEVTOOLS);
    } catch (e) {
      console.error("Vue Devtools failed to install:", e.toString());
    }
  }
  createWindow();
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
