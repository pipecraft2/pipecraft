module.exports = {
  moduleFileExtensions: ["js", "json", "vue"],
  transform: {
    "^.+\\.js$": "babel-jest",
    "^.+\\.vue$": "vue-jest",
  },
  testMatch: ["**/tests/unit/**/*.spec.js"],
  testPathIgnorePatterns: ["<rootDir>/dist_electron/"],
  globals: {
    "vue-jest": {
      resources: {
        scss: ["./src/assets/styles/_variables.scss"],
      },
    },
  },
  moduleNameMapper: {
    "^vue$": "vue/dist/vue.common.js",
    "^@/(.*)$": "<rootDir>/src/$1",
    "^envfile$": require.resolve("envfile"),
    "\\.(css|less|scss|sass)$": "identity-obj-proxy", // Mock CSS imports
  },
  testEnvironment: "jsdom", // Use jsdom environment
  setupFiles: ["<rootDir>/tests/setup.js"],
};
