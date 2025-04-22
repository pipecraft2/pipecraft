# PipeCraft2 <img src='src/assets/PipeCraft2_logo.png' align="right" height="180" />

PipeCraft2 is a user-friendly GUI software for metabarcoding data analysis. It provides:
- Ready-to-run pipelines for common metabarcoding workflows
- Individual analysis modules for custom workflows
- Integration with popular bioinformatics tools

## Quick Start

### For Users

📥 [Download Latest Release (v1.1.0)](https://github.com/pipecraft2/pipecraft/releases/tag/v1.1.0)

📚 [User Guide](https://pipecraft2-manual.readthedocs.io/en/latest/)

### For Developers

Pre-built Docker images available on [DockerHub](https://hub.docker.com/u/pipecraft).

📚 [Developer Guide](https://pipecraft2-manual.readthedocs.io/en/1.0.0/for_developers.html)

#### Prerequisites

- [NodeJS](https://nodejs.org/en/download/) (make sure to install `build tools` on Windows, or `build-essential` on Ubuntu)
- [Yarn package manager](https://classic.yarnpkg.com/en/docs/install/#windows-stable)
- [Docker](https://www.docker.com/get-started)
- [Git](https://git-scm.com/downloads)

#### Installation

```bash
# Clone the repository
git clone https://github.com/pipecraft2/pipecraft
cd pipecraft

# Install dependencies and setup PipeCraft
yarn run install_pipe

# Start PipeCraft in development mode
yarn electron:serve
```

Linux:


### Install Build Dependencies
These packages are required for compiling native modules and building the application:

```bash
sudo apt update
sudo apt install python3-dev python3-pip python3-setuptools
sudo apt install build-essential
```
# Download and install nvm:
```bash
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.2/install.sh | bash
```
# Restart the shell
```bash
\. "$HOME/.nvm/nvm.sh"
```
# Download and install Node.js:
```bash
nvm install 16
```
```bash
npm install --global yarn
```
```bash
git clone https://github.com/pipecraft2/pipecraft
```
```bash
cd pipecraft
```
# Install dependencies
```bash
yarn install
```