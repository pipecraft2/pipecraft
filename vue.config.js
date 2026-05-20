module.exports = {
  transpileDependencies: ["vuetify"],
  pluginOptions: {
    electronBuilder: {
      nodeIntegration: true,
      externals: ["dockerode", "node-pty-prebuilt-multiarch"],
      builderOptions: {
        // Explicit GitHub target so CI publishes without relying on git remote metadata.
        publish: [
          {
            provider: "github",
            owner: "pipecraft2",
            repo: "pipecraft",
            // Must match GitHub: if a published release already exists for the tag,
            // electron-builder cannot publish as draft (see GH "existing type not compatible").
            releaseType: "release",
          },
        ],
        // Windows: NSIS enables electron-updater; portable builds do not support auto-update.
        win: {
          icon: "build/icon.ico",
          target: [
            {
              target: "nsis",
              arch: ["x64"],
            },
          ],
        },
        // Linux: single AppImage for releases. FUSE-free option: APPIMAGE_EXTRACT_AND_RUN=1.
        linux: {
          target: ["AppImage"],
          icon: "build/icons",
          category: "Science",
          maintainer: "PipeCraft2 Team: martin.metsoja@ut.ee, sten.anslan@ut.ee",
          vendor: "PipeCraft2",
          synopsis:
            "Bioinformatics application that implements various popular tools for metabarcoding data analyses.",
          description: "PipeCraft is a desktop application for metabarcoding data analysis.",
          desktop: {
            Name: "PipeCraft2",
            Comment: "Software for metabarcoding data analysis",
            Categories: "Science;Biology;Bioinformatics",
            Terminal: false,
            Type: "Application",
            StartupNotify: true,
            StartupWMClass: "pipecraft",
          },
        },
        appImage: {
          artifactName: "${productName}-${version}-linux-${arch}.AppImage",
        },
        mac: {
          target: [
            {
              target: "dmg",
              // One universal binary + single DMG: Intel (x64) and Apple Silicon (arm64)
              // without building two separate DMGs in one job (avoids intermittent
              // GitHub runner hdiutil "Resource busy" when producing arm64 + x64 images).
              arch: ["universal"],
            },
          ],
          icon: "build/icon.icns", 
          hardenedRuntime: true,      // Required for macOS 10.15+ (Catalina)
          gatekeeperAssess: false,    // Skip Gatekeeper assessment
          entitlements: "build/entitlements.mac.plist",        // Path to entitlements
          entitlementsInherit: "build/entitlements.mac.plist", // Child process entitlements
        },
        afterSign: "build/notarize.js",
        appx: {
          applicationId: "pipecraft",
        },
        appId: "pipecraft",
        productName: "pipecraft",
        extraResources: [
          "src/pipecraft-core/**/*",
          "!src/pipecraft-core/**/.git/**",
          "!src/pipecraft-core/**/.gitignore",
          "!src/pipecraft-core/**/.gitattributes"
        ],
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
