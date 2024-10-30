module.exports = {
  moduleFileExtensions: ["js", "json", "vue"],
  transform: {
    "^.+\\.js$": "babel-jest",
    "^.+\\.vue$": "vue-jest",
  },
  testMatch: ["**/tests/unit/**/*.spec.js"],
  setupFiles: ["<rootDir>/tests/setup.js"],
};
