import Vue from "vue";
import Vuetify from "vuetify"; // Ensure it's imported

Vue.use(Vuetify);

// Mock electron's ipcRenderer

global.electron = {
  ipcRenderer: {
    send: jest.fn(),
    on: jest.fn(),
  },
};

// Mock @electron/remote for Jest
jest.mock("@electron/remote", () => ({
  dialog: {
    showOpenDialog: jest.fn(),
  },
  require: jest.fn(() => ({
    spawn: jest.fn(() => ({
      on: jest.fn(),
    })),
  })),
}));

// Mock Dockerode
jest.mock("dockerode", () => {
  return jest.fn().mockImplementation(() => {
    return {
      run: jest.fn().mockResolvedValue([{ StatusCode: 0 }, {}]),
    };
  });
});
