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
          maintainer: "PipeCraft2 Team: martin.metsoja@ut.ee, sten.anslan@ut.ee",
          vendor: "PipeCraft2",
          synopsis: "Bioinformatics application that implements various popular tools for metabarcoding data analyses.",
          description: "PipeCraft is a desktop application for metabarcoding data analysis.",
          desktop: {
            Name: "PipeCraft2",
            Comment: "Software for metabarcoding data analysis",
            Categories: "Science;Biology;Bioinformatics",
            Terminal: false,
            Type: "Application",
            StartupNotify: true,
            StartupWMClass: "pipecraft"
          }
        },
        mac: { 
          target: [
            {
              target: "dmg",
              arch: ["arm64", "x64" ]  // Architectures specified in target
            }
          ],
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
        files: [
          "**/*",
          "!**/node_modules/*/{CHANGELOG.md,README.md,README,readme.md,readme}",
          "!**/node_modules/*/{test,__tests__,tests,powered-test,example,examples}",
          "!**/node_modules/*.d.ts",
          "!**/node_modules/.bin",
          "!**/*.{iml,o,hprof,orig,pyc,pyo,rbc,swp,csproj,sln,xproj}",
          "!.editorconfig",
          "!**/._*",
          "!**/{.DS_Store,.git,.hg,.svn,CVS,RCS,SCCS,.gitignore,.gitattributes}",
          "!**/{__pycache__,thumbs.db,.flowconfig,.idea,.vs,.nyc_output}",
          "!**/{appveyor.yml,.travis.yml,circle.yml}",
          "!**/{npm-debug.log,yarn.lock,.yarn-integrity,.yarn-metadata.json}",
          {
            from: "node_modules/node-pty-prebuilt-multiarch/build/Release/",
            to: "node_modules/node-pty-prebuilt-multiarch/build/Release/",
            filter: ["**/*"]
          }
        ]
      },
    },
  },
};
