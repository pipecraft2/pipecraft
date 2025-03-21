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
          target: [
            {
              target: "portable",
              arch: ["x64"]
            }
          ]
        },
        linux: {
          target: ["deb", "rpm", "AppImage"],  // Build both DEB and RPM packages
          icon: "build/icons",     // Directory containing multiple icon sizes
          category: "Science",
          maintainer: "PipeCraft Team martin.metsoja@ut.ee",
          vendor: "PipeCraft",
          synopsis: "Metabarcoding application",
          description: "PipeCraft is a desktop application for metabarcoding data analysis.",
          desktop: {
            Name: "Pipecraft",
            Comment: "Metabarcoding application",
            Categories: "Science;Biology;",
            Terminal: true,
            Type: "Application",
            StartupNotify: true,
            StartupWMClass: "pipecraft"
          }
        },
        mac: { 
          target: "default", 
          icon: "build/icon.icns", 
          hardenedRuntime: true,      // Required for macOS 10.15+ (Catalina)
          gatekeeperAssess: false,    // Skip Gatekeeper assessment
          entitlements: "build/entitlements.mac.plist",        // Path to entitlements
          entitlementsInherit: "build/entitlements.mac.plist", // Child process entitlements
          identity: "MARTIN METSOJA (3CYD3SH29Q)", // Your signing identity 
        },
        afterSign: "build/notarize.js",
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
