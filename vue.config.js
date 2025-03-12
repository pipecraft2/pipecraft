module.exports = {
  transpileDependencies: ["vuetify"],
  pluginOptions: {
    electronBuilder: {
      nodeIntegration: true,
      externals: ["dockerode", "node-pty-prebuilt-multiarch"],
      builderOptions: {
        publish: ["github"],
        win: {
          icon: "build/icon.ico",
          // target: [
          //   {
          //     target: "portable",
          //     arch: ["x64"]
          //   }
          // ]
        },
        linux: {
          target: "AppImage",
          icon: "build/icon.png",
        },
        mac: { target: "default", icon: "build/icon.icns" },
        appx: {
          applicationId: "pipecraft",
        },
        appId: "pipecraft",
        productName: "pipecraft",
        extraResources: ["src/pipecraft-core"],
      },
    },
  },
};
